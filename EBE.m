, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Biomass Estimation EBE - version 2.0
%
% Biomass mapping and budgeting module
%
% Melissa Latella & Tommaso Raimondo
%
% Created: 01 November 2021
% Last update: 26 July 2022
%
% melissa.latella@polito.it
%
% Cite as: "Latella, M., Raimondo, T., Belcore, L., Salerno, L.,
%                       and Camporeale, C. (2022), On the integration of
%                       LiDAR and field data for riparian biomass
%                       estimation, Journal of Environmental Management".
%
%Aim: Mapping riparian vegetation biomass considering trees and shrubs along river
%        channel(s) and floodplain by combining field data and LiDAR point clouds.
%
% Method: The algorithm defines the vegetated areas according to input RGB 
%           images and shapefiles indicating the water mask and the Region 
%           of Interest (RoI). 
%         The algorithm maps the biomass of the trees standing within
%           the RoI by using two different methods: (i) the trees outside 
%           the Shrub Area (SA) or on islands are segmented according to 
%           the watershed algorithm by Chen et al. (2006), using as 
%           reference treetops the individual trees previously detected by 
%           the EBE_ITDM.m module (Latella et al., 2021). (ii) within the
%           SA, the tree crowns are segmented according to Hasenauer
%           (1997) crown width. The biomass is always computed according to
%           field-based allometric relationships (Latella et al., 2020, 
%           Latella et al., 2021)
%         The algorithm maps the shrub biomass by detecting individuals or
%           clusters of shrubs based on the mask previously generated.
%           (The tree crowns are removed from the mask). For each detected
%           shrub, it clips the point clouds and computes its volume by 
%           following the concave hull slicing method proposed by
%           Yan et al. (2019) and Moreira and Santos (2007). It converts
%           volumes into biomass through a bulk density coefficients
%           estimated from field measurements.
%          The output consists of two json shapefiles for each input area,
%           one concerning plant estimated biomass and another for the
%           associated error. The json files can be visualized and
%           processed in GIS environment to obtain maps and budgets.
%          Optionally, an overall json file might be generated if more than
%           one area is analyzed.
%
% Dependencies:
% 0. Matlab version: 2021 or followings.
% 1. output from EBE_ITDM.m by Melissa Latella.
% 2. EBE_concave_hull.m by Tommaso Raimondo & Melissa Latella. In turn,
% concave_hull demands the function selfintersect.m by Antoni J. Canós 
% (2022). Fast and Robust Self-Intersections 
% (https://www.mathworks.com/matlabcentral/fileexchange/13351-fast-and-robust-self-intersections), 
% MATLAB Central File Exchange. Retrieved July 26, 2022.
% 3. jsonInitialize.m by Melissa Latella.
% 4. jsonUpdate.m by Melissa Latella
% 5. Lidar Toolbox™ by Matlab.
% 6. Parallel Toolbox™ by Matlab (otherwise substitute "parfor" with "for"
%    in the script and comment rows from 156 to 160, number of cores row 84).
%
% Warning (i): rows 140 to 152 and 283 might need to be modified depending 
% on the format of the input filenames.
%
% Warning (ii): rows 756 to 781 are commented. They perform an optional 
% buffer of shrub clusters and then merge the overlapping buffer. 
% It is time consuming.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
warning('off','all')


%%%%%%%%%%%%%%% parameter setup %%%%%%%%%%%%%%%

%options
existingClippedClouds = false;
writeTotalJson = false;
numberCores = 8;

%RGB parameters
redLowerThreshold = 0;
redUpperThreshold = 0.25;
greenLowerThreshold = 0;
greenUpperThreshold = 0.15;
blueLowerThreshold = 0;
blueUpperThreshold = 0.1;

%height threshold to discriminate shrubs from trees
treeHeightThreshold = 2.5;%m, based on field observations

%height threshold to discriminate ground from vegetation
groundThreshold = 0.10; %m, based on field observations

%pixel threshold of minimum cluster size
clusterThreshold = 4; %pixels (in the example the cell side is 0.33 m)

%LiDAR point density
lidarPointsDensity = 9; %points/m2, from acquisition specifications
tileColumnsStartFrom = 'south';
tileRowsStartFrom =  'west';

%buffer radius
clusterRadius = 1; %m

%concave hull function parameter
dh = 0.30; %m, initial class height
pointThreshold = 5; %minimum point per class

%vegetation bulk density
bulkDensity = 2.578; %kg/m3, weighted average density based on field data

%volume correction factor
volumeCorrection = 2.25; %-, based on field data

%allometric coefficients for different tree species
fullName = {'Populus spp.';'Quercus spp.';'Robinia Pseudoacacia'};
shortName = {'populus';'quercus';'robinia'};
aParameter = [1.2560; 4.8273; 3.6539];%based on field data
bParameter = [0.7211; 0.3772; 0.4496];%based on field data
greenWoodDensity = [0.78; 1.05; 1.01].*10^3;%kg/m3, from the literature
species = table(fullName, shortName, aParameter, bParameter, greenWoodDensity);

%error statistics
treeDensityStd = 200; %kg/m3
AStd = [0.15841; 0.62769; 1.1935];
BStd = [0.036968; 0.038863; 0.11088];
HStd = 1.21; %m
shrubDensityStd = 1.1336; %kg/m3
ALSVolStd = 0.1673; %m3


%%%%%%%%%%%%%%% input and output definition %%%%%%%%%%%%%%%

%json properties
jsonName = 'CNR2021';
jsonCrs = '32632';

%input files
inputFolder = 'Input';
outputFolder = 'Output';
shrubPolygonFileName = 'ShrubArea.poly';%text format e.g., .poly .txt
waterPolygonFileName = 'WaterMask.shp';%.shp
islandFileName = 'Islands.poly';%text format e.g., .poly .txt
roiPolygonFileName = 'RegionOfInterest.poly';%text format e.g., .poly .txt
treeDirectory = strcat(inputFolder,'/ITDM_Results/');
jsonFileName = 'EBE';


%%%%%%%%%%%%%%% parallel definition %%%%%%%%%%%%%%%

delete(gcp('nocreate'))
poolobj = parpool('local',numberCores);
spmd
    warning('off','MATLAB:rmpath:DirNotFound')
end


%%%%%%%%%%%%%%% directory creation %%%%%%%%%%%%%%%

%message
fprintf('%s\n', 'START')

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
end
if ~exist(strcat(outputFolder, '/Trees'), 'dir')
    mkdir(strcat(outputFolder, '/Trees'))
end
if ~exist(strcat(outputFolder, '/Shrubs'), 'dir')
    mkdir(strcat(outputFolder, '/Shrubs'))
end
if ~exist(strcat(outputFolder,'/exportedJson'), 'dir')
    mkdir(strcat(outputFolder,'/exportedJson'))
end
if ~exist(strcat(outputFolder,'/exportedJson/Single'), 'dir')%results for each input area
    mkdir(strcat(outputFolder,'/exportedJson/Single'))
end

%%%%%%%%%%%%%%% input retrieval %%%%%%%%%%%%%%%

%%%%input

%forest types
forestList = struct2cell(dir(strcat(inputFolder,'/ForestTypePolygons/*.txt')));
forestNumber = size(dir(strcat(inputFolder,'/ForestTypePolygons/*.txt')),1);
for i = 1 : forestNumber
    speciesName = cell2mat(forestList(1,i));
    polygonFileName = strcat(inputFolder,'/ForestTypePolygons/',speciesName);
    polygonRead = readmatrix(polygonFileName, 'FileType', 'text');
    eval(['forestPolygon' num2str(i) '= polygonRead;']);
    eval(['forestPolyshape' num2str(i) '= polyshape(forestPolygon' num2str(i) '(:,1), forestPolygon' num2str(i) '(:,2));']);
end

%RGB images and the associated point clouds
rgbList = struct2cell(dir(strcat(inputFolder,'/RGBs/*.tif')));
pointcloudList = struct2cell(dir(strcat(inputFolder,'/PointClouds/*.las')));
rgbNames = sort(rgbList(1,:));
pointcloudsNames = sort(pointcloudList(1,:));
%number of input areas
areaNumber = size(rgbList,2);

%region of interest (RoI)
roiPolygon0 = readmatrix(strcat(inputFolder,'/Polygons/',roiPolygonFileName), 'FileType', 'text');
roiPolygonX = roiPolygon0(:,1);
roiPolygonY = roiPolygon0(:,2);
roiPolygon = polyshape(roiPolygonX,roiPolygonY);

%shrub area (SA)
saPolygon0 = readmatrix(strcat(inputFolder,'/Polygons/',shrubPolygonFileName), 'FileType', 'text');
polygonX = saPolygon0(:,1);
polygonY = saPolygon0(:,2);
saPolygon = polyshape(polygonX,polygonY);

%water mask
w = shaperead(strcat(inputFolder,'/Polygons/',waterPolygonFileName));
waterPolygon = [w.X' w.Y'];

%islands
islandPolygon = readmatrix(strcat(inputFolder,'/Polygons/',islandFileName), 'FileType', 'text');

%%%%first operations

%separate wet regions
waterNumber = numel(find(isnan(waterPolygon(:,1))))+1;
waterPolyshape = polyshape();
inds = find(isnan(waterPolygon(:,1)));
indEnd = 1;
for i = 1 : waterNumber
    indStart = indEnd;
    if i == waterNumber
        indEnd = numel(waterPolygon(:,1));
    else
        indEnd = inds(i)-1;
    end
    waterPolyshape(i,:) = polyshape( waterPolygon(indStart:indEnd,1)', waterPolygon(indStart:indEnd,2)');
end

%separate islands
islandNumber = numel(find(isnan(islandPolygon(:,1))))+1;
islandPolyshape = polyshape();
inds = find(isnan(islandPolygon(:,1)));
for i = 1 : islandNumber
    if i == 1
        indStart = 1;
    end
    if i == islandNumber
        indEnd = numel(islandPolygon(:,1));
    else
        indEnd = inds(i);
    end
    islandPolyshape(i,:) = polyshape( islandPolygon(indStart:indEnd,1)', islandPolygon(indStart:indEnd,2)');
    indStart = indEnd;
end


%%%%%%%%%%%%%%% biomass computation %%%%%%%%%%%%%%%

%%%%output initialization
if writeTotalJson
    jsonBiomTot = jsonInitialize(jsonName,jsonCrs);
    jsonErrTot = jsonInitialize(jsonName,jsonCrs);
end

%%%%loop over the file list

for s = 1 : areaNumber

    %message
    fprintf('Tile: %s\n', cell2mat(rgbNames(1,s)))

    %%%%file import
    rgbFileName = cell2mat(rgbNames(1,s));
    [~,tileName,~] = fileparts(rgbFileName);
    rgbImage = imread((strcat(inputFolder,'/RGBs/',rgbFileName)));
    treeFileName = dir(strcat(treeDirectory,tileName, '*.txt'));
    treeFileName = treeFileName.name;
    trees = readmatrix(strcat(treeDirectory,treeFileName));

    %message
    fprintf('%s\n', 'RGB Input: defined')


    %%%%adjust histogram of rgb images
    rgbNorm = double(rgbImage).^4;
    rgbNorm = rgbNorm/max(rgbNorm(:)); %now RGB within [0,1]


    %%%%get coordinates from the RGB image
    [Z,R] = readgeoraster ((strcat(inputFolder,'/RGBs/',rgbFileName)));
    info = georasterinfo((strcat(inputFolder,'/RGBs/',rgbFileName)));
    missingDataInd = info.MissingDataIndicator;
    Z = standardizeMissing(Z,missingDataInd);
    [rowNumber,colNumber] = size(Z,1:2);
    [rows,cols] = ndgrid(1:rowNumber,1:colNumber);
    [tileX,tileY] = intrinsicToWorld(R,cols,rows);


    %%%%define a mask of the RoI

    %bounding box of the tile
    bboxX = [tileX(1,1,1); tileX(1,size(tileX,2),1); tileX(size(tileX,1),size(tileX,2),1); tileX(size(tileX,1),1,1);tileX(1,1,1)];
    bboxY = [tileY(1,1,1); tileY(1,size(tileY,2),1); tileY(size(tileY,1),size(tileY,2),1); tileY(size(tileY,1),1,1);tileY(1,1,1)];
    bbox = polyshape(bboxX,bboxY);

    %intersection between bounding the RoI and the tile
    RoIntersect = intersect(roiPolygon, bbox);

    %if the RoI overlaps the RGB image
    if ~isempty(RoIntersect.Vertices)

        %initialize json
        clear jsonBiom jsonErr
        [~,fileName,~] = fileparts(rgbFileName);
        jsonName = fileName;
        jsonBiom = jsonInitialize(jsonName,jsonCrs);
        jsonErr = jsonInitialize(jsonName,jsonCrs);

        %roi mask
        roiMask = roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , RoIntersect.Vertices(~isnan(RoIntersect.Vertices(:,1)),1), RoIntersect.Vertices(~isnan(RoIntersect.Vertices(:,1)),2));%roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , polygonX, polygonY);

        %%%%define a mask of the vegetated areas based on RGB properties

        %separate RGB bands
        redBand = rgbNorm (:,:,1);
        greenBand = rgbNorm (:,:,2);
        blueBand = rgbNorm (:,:,3);

        %create a mask per band
        redMask = (redBand > redLowerThreshold & redBand < redUpperThreshold);
        greenMask = (greenBand > greenLowerThreshold & greenBand < greenUpperThreshold);
        blueMask = (blueBand > blueLowerThreshold & blueBand < blueUpperThreshold);

        %%%%define a mask of the vegetated areas based on RGB properties

        %comprehensive RGB mask
        rgbMask = uint8(redMask & greenMask & blueMask);
        vegetationMask = (rgbMask & roiMask);

        %%%%define a mask of the active channel
        for i = 1 : waterNumber
            waterIntersect =  intersect(waterPolyshape(i,:), bbox);

            %%if the tile overlaps the wet area
            if ~isempty(waterIntersect.Vertices)
                intPoly = intersect(bbox,waterPolyshape(i,:));

                %check polygon division
                inds = find(isnan(intPoly.Vertices(:,1)));
                indEnd = 1;
                if ~isempty(inds)
                    for j = 1:numel(inds)+1
                        indStart = indEnd;
                        if j == numel(inds)+1
                            indEnd = numel(intPoly.Vertices(:,1));
                        else
                            indEnd = inds(j)-1;
                        end
                        intPolyTemp = polyshape( intPoly.Vertices(indStart:indEnd,1)', intPoly.Vertices(indStart:indEnd,2)');
                        waterPolygonX =  intPolyTemp.Vertices(~isnan(intPolyTemp.Vertices(:,1)),1);
                        waterPolygonY =  intPolyTemp.Vertices(~isnan(intPolyTemp.Vertices(:,1)),2);
                        waterMask = roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , waterPolygonX, waterPolygonY); %don't change the order: it may rotate the mask

                        %update the vegetation mask
                        vegetationMask = vegetationMask - waterMask;
                        vegetationMask(vegetationMask==-1) = 0;

                    end

                else
                    waterPolygonX = intPoly.Vertices(~isnan(intPoly.Vertices(:,1)),1);
                    waterPolygonY = intPoly.Vertices(~isnan(intPoly.Vertices(:,1)),2);
                    waterMask = roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , waterPolygonX, waterPolygonY); %don't change the order: it may rotate the mask

                    %update the vegetation mask
                    vegetationMask = vegetationMask - waterMask;
                    vegetationMask(vegetationMask==-1) = 0;

                end

            end%if tile overlaps the wet area

        end%loop over wet areas

        %%intersection between bounding box vegetated islands
        islandInside = [];
        for i = 1 : islandNumber
            islandIntersect =  intersect(islandPolyshape(i,:), bbox);
            if ~isempty(islandIntersect.Vertices)
                islandPolygonX = islandPolyshape(i,:).Vertices(~isnan(islandPolyshape(i,:).Vertices(:,1)),1);
                islandPolygonY = islandPolyshape(i,:).Vertices(~isnan(islandPolyshape(i,:).Vertices(:,1)),2);

                %water mask
                islandMask = roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , islandPolygonX, islandPolygonY); %don't change the order: it may rotate the mask!

                %update the vegetation mask
                vegetationMask = vegetationMask - islandMask;
                vegetationMask(vegetationMask==-1) = 0;

                %polygon index inside the tile
                islandInside = [islandInside,i];

            end
        end

        %message
        fprintf('%s\n', 'Roi, water and vegetation masks: created')

        %%trees inside the SA

        %if the tile intersects the shrub area
        saIntersect = intersect(saPolygon, bbox);
        if ~isempty(saIntersect.Vertices)
            saMask = roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , polygonX, polygonY); %don't change the order: it may rotate the mask!

            %update the vegetation mask
            vegetationMask = vegetationMask.* saMask;
            vegetationMask(vegetationMask==-1) = 0;

            %find trees in the shrub area
            treeInIndex = inpolygon(trees(:,1), trees(:,2), saPolygon0(:,1), saPolygon0(:,2));
            treeIn = trees(treeInIndex == 1,:);

            %remove trees on the islands
            treeIslandsIndex = [];
            if ~isempty(islandInside)
                for h = islandInside(1) : islandInside(end)
                    treeIslandsIndex = [treeIslandsIndex;find(inpolygon(treeIn(:,1),treeIn(:,2),islandPolyshape(h,:).Vertices(:,1),islandPolyshape(h,:).Vertices(:,2)))];
                end
                treeIn(treeIslandsIndex,:) = [];
            end

            %detect trees higher than treeHeightThreshold
            treeOverIndex = find(treeIn(:,3) > treeHeightThreshold);

            %save detected tree list (i.e., trees above threshold in the SA)
            treeOutput = treeIn(treeOverIndex,:);
            treeOutput = array2table(treeOutput);
            treeOutput.Properties.VariableNames = {'X','Y','Z'};
            [~,fileName,~] = fileparts(rgbFileName);
            treeOutputFileName = strcat(outputFolder,'/Trees/', 'T_', fileName,'.txt');
            writetable(treeOutput, treeOutputFileName);

            %%remove the points in the cloud corresponding to the excluded trees
            treePolygon = polyshape();
            treeOver = treeIn(treeOverIndex,:);

            %segmentate tree in the SA
            for i = 1:size(treeOver ,1)

                %tree's crown feature through Hasenauer Equation
                %Hasenauer, H. (1997). Dimensional relationships of open-grown trees in Austria.
                %Forest Ecology and Management, 96(3), 197-206.
                crownWidth = exp(-0.3526 + 0.8680*log(treeOver (i,3))); %m  Hasenauer equation for Salicaceus species
                treeRadius = crownWidth/2; % m
                crownArea = pi*treeRadius^2;

                %define buffer around the tree with a radius == treeRadius
                treePolygon(i,:) = polybuffer(treeOver (i,1:2),'line',treeRadius);

                %remove the buffered areas from the vegetation mask
                treeMask = roipoly([min(bboxX); max(bboxX)], [max(bboxY); min(bboxY)], rgbNorm , treePolygon(i,:).Vertices(:,1), treePolygon(i,:).Vertices(:,2));%dont' change the order
                vegetationMask = vegetationMask - treeMask;
                vegetationMask(vegetationMask==-1) = 0;

                %find tree species
                forestTypeInd = 0;
                for j = 1 : forestNumber
                    eval(['typeIndex = inpolygon( treeOver(i,1), treeOver(i,2), forestPolygon' num2str(j) '(:,1), forestPolygon' num2str(j) '(:,2));']);
                    if numel(typeIndex)  > 0
                        forestTypeInd = j;
                        break
                    end
                end%forest type
                if forestTypeInd == 0
                    forestTypeInd = 1; %generic = populus
                end

                %estimate biomass of individual trees through species-dependent allometric relationships
                treeHeight  = treeOver (i,3); %m
                treeDiameter  = 0.01 * (treeHeight/aParameter(forestTypeInd))^(1/bParameter(forestTypeInd)); %m
                treeBiomass  = (0.25 * pi * treeDiameter^2) * treeHeight * greenWoodDensity(forestTypeInd); %kg

                %estimate the associated error
                treeVariance = (treeBiomass^2)*( (treeDensityStd/greenWoodDensity(forestTypeInd))^2 + (-2*AStd(forestTypeInd)/(aParameter(forestTypeInd)*bParameter(forestTypeInd)))^2 + (-2*log(treeHeight/aParameter(forestTypeInd))*BStd(forestTypeInd)/(bParameter(forestTypeInd)^2))^2 + ((bParameter(forestTypeInd)+2)*HStd/(bParameter(forestTypeInd)*treeHeight))^2 );
                treeStd = sqrt(treeVariance);

                %write json file
                listingtrees = cell({i});
                filenameTrees = strcat('T_',string(cell2mat(listingtrees)));
                jsonBiom = jsonUpdate(jsonBiom, 'tree', filenameTrees, treeBiomass,treePolygon(i,:));
                jsonErr = jsonUpdate(jsonErr, 'tree', filenameTrees, treeStd,treePolygon(i,:));
                if writeTotalJson
                    jsonBiomTot = jsonUpdate(jsonBiomTot, 'tree', filenameTrees, treeBiomass,treePolygon(i,:));
                    jsonErrTot = jsonUpdate(jsonErrTot, 'tree',filenameTrees, treeStd,treePolygon(i,:));
                end

            end%tree within the SA

        end%if the SA intersects the tile


        %%trees outside the SA but inside RoI

        %read point cloud file and related information
        ptCloudFileName = cell2mat(pointcloudsNames(1,s));
        lasReader = lasFileReader(strcat(inputFolder,'/PointClouds/', ptCloudFileName));
        [ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","ScanAngle");
        groundIdx = segmentGroundSMRF(ptCloud);
        groundPtCloud = select(ptCloud,groundIdx);
        groundPoints = groundPtCloud.Location;

        %remove duplicate points along x- and y-axes
        [zUnique,xyUnique] = groupsummary(groundPoints(:,3),groundPoints(:,1:2),@mean);
        xyUnique = [xyUnique{:}];

        %normalize point elevation
        interpFunction = scatteredInterpolant(double(xyUnique),double(zUnique),"natural");
        ptElevation = interpFunction(double(ptCloud.Location(:,1)),double(ptCloud.Location(:,2)));
        ptCloudNormalized = ptCloud.Location;
        ptCloudNormalized(:,3) = ptCloudNormalized(:,3) - ptElevation;

        %generate the canopy height model (CHM)
        gridResolution = 1/sqrt(lidarPointsDensity); %optimal grid resolution according to the literature
        canopyModel = pc2dem(pointCloud(ptCloudNormalized),gridResolution,"CornerFillMethod","max");
        canopyModel(isnan(canopyModel) | canopyModel<0) = 0;
        chmFilter = fspecial("gaussian",[5 5],1);
        canopyModel = imfilter(canopyModel,chmFilter,'replicate','same');

        %trees in intrinsic coordinates
        treeXWorld = trees(:,1);
        treeYWorld = trees(:,2);
        R.CellExtentInWorldX = gridResolution;
        R.CellExtentInWorldY = gridResolution;
        R.ColumnsStartFrom = tileColumnsStartFrom;
        R.RowsStartFrom =  tileRowsStartFrom;
        treeYIntrinsic = round(((treeYWorld - ptCloud.YLimits(1))/gridResolution) + 1);%rows
        treeXIntrinsic = round(((treeXWorld - ptCloud.XLimits(1))/gridResolution) + 1);%columns

        %transfer labels to each point
        rowId = ceil((ptCloud.Location(:,2) - ptCloud.YLimits(1))/gridResolution) + 1;
        colId = ceil((ptCloud.Location(:,1) - ptCloud.XLimits(1))/gridResolution) + 1;       
        minTreeHeight = treeHeightThreshold;
        markerImage = false(size(canopyModel));
        validTreeTops = sub2ind(size(canopyModel),treeYIntrinsic,treeXIntrinsic);
        markerImage(validTreeTops) = true;
        canopyComplement = -canopyModel;
        minImage = imimposemin(canopyComplement, markerImage);
        label2D = watershed(minImage, 8);
        label2D(canopyModel < minTreeHeight) = 0;
        label2D = bwlabel(label2D ~= 0, 8);
        ind = sub2ind(size(label2D),rowId,colId);
        label3D = label2D(ind);

        %extract valid labels and corresponding points
        validSegIds = label3D ~= 0;
        veglabel3D = label3D(validSegIds);

        %remove tree outside
        %segmentate only inside RoI but outside SA
        treeNumber = max(veglabel3D);
        ptCloud = double(ptCloud.Location(:,:));

        %find trees inside on the islands
        if ~isempty(islandInside)

            for j = 1:treeNumber

                for h = islandInside(1) : islandInside(end)

                    %find trees (i.e., cells of the rasterized area) in the RoI, on the islands and in the SA
                    cellIndex = find(label3D == j);
                    cellRoiInd = inpolygon(ptCloudNormalized(cellIndex,1),ptCloudNormalized(cellIndex,2),roiPolygonX,roiPolygonY);
                    cellIslandsInd = inpolygon(ptCloudNormalized(cellIndex,1),ptCloudNormalized(cellIndex,2),islandPolyshape(h,:).Vertices(:,1),islandPolyshape(h,:).Vertices(:,2));
                    cellSaInd = inpolygon(ptCloudNormalized(cellIndex,1),ptCloudNormalized(cellIndex,2),polygonX,polygonY);

                    %find trees on islands OR in the RoI but outside the SA
                    if sum(cellIslandsInd) ~= 0 || (sum(cellRoiInd) ~= 0 && sum(cellSaInd) == 0) 
                        ptCloudRoi = double(ptCloudNormalized(cellIndex,1:2));
                        [treeBoundaryVertices,~] = boundary(ptCloudRoi(:,1:2),0.5);

                        %find tree species (i.e., if points across
                        %different species, choose where most points are)
                        forestTypeInd = 0;
                        forestCount = zeros(forestNumber,1);
                        for v = 1 : forestNumber
                          treePolygonTemp = polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2));
                          eval(['forestCount(v,1) = area(intersect( treePolygonTemp, forestPolyshape' num2str(v) '));']);
                          [~, forestTypeInd ] = max(forestCount);%it choose the first by default if more maxs
                        end%forest type
                        if sum(forestCount) == 0 %if points do not belonging to a known forest type
                            forestTypeInd = 1; %generic = populus
                        end

                        %estimate biomass of individual trees through species-dependent
                        %allometric relationships
                        treeHeight  = max(ptCloud(cellIndex,3)) - min(ptCloud(cellIndex,3)); %m
                        treeDiameter  = 0.01 * (treeHeight/aParameter(forestTypeInd))^(1/bParameter(forestTypeInd)); %m
                        treeBiomass  = (0.25 * pi * treeDiameter^2) * treeHeight * greenWoodDensity(forestTypeInd); %kg

                        %estimate the associated error
                        treeVariance = (treeBiomass^2)*( (treeDensityStd/greenWoodDensity(forestTypeInd))^2 + (-2*AStd(forestTypeInd)/(aParameter(forestTypeInd)*bParameter(forestTypeInd)))^2 + (-2*log(treeHeight/aParameter(forestTypeInd))*BStd(forestTypeInd)/(bParameter(forestTypeInd)^2))^2 + ((bParameter(forestTypeInd)+2)*HStd/(bParameter(forestTypeInd)*treeHeight))^2 );
                        treeStd = sqrt(treeVariance);

                        %write json file
                        listingtrees = cell({j});
                        filenameTrees = strcat('T_',string(cell2mat(listingtrees)));
                        jsonBiom = jsonUpdate(jsonBiom, 'tree', filenameTrees, treeBiomass,polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                        jsonErr = jsonUpdate(jsonErr, 'tree', filenameTrees, treeStd,polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                        if writeTotalJson
                            jsonErrTot = jsonUpdate(jsonBiomTot, 'tree', filenameTrees, treeBiomass, polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                            jsonErrTot = jsonUpdate(jsonErrTot, 'tree', filenameTrees, treeStd, polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                        end

                    end

                end %loop over islands

                ptCloud(cellIndex,1:2) = NaN;

            end %loop over trees

        else

            for j = 1:treeNumber

                %find trees in the RoI and in the SA
                cellIndex = find(label3D == j);
                cellRoiInd = inpolygon(ptCloudNormalized(cellIndex,1),ptCloudNormalized(cellIndex,2),roiPolygonX,roiPolygonY);
                cellSaInd = inpolygon(ptCloudNormalized(cellIndex,1),ptCloudNormalized(cellIndex,2),polygonX,polygonY);

                %find trees in the RoI but outside the SA
                if sum(cellRoiInd) ~= 0 && sum(cellSaInd) == 0
                    ptCloudRoi = double(ptCloudNormalized(cellIndex,1:2));
                    [treeBoundaryVertices, ~] = boundary(ptCloudRoi(:,1:2),0.5);

                    %find tree species (i.e., if points across
                    %different species, choose where most points are)
                    forestTypeInd = 0;
                    forestCount = zeros(forestNumber,1);
                    for v = 1 : forestNumber
                        treePolygonTemp = polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2));
                        eval(['forestCount(v,1) = area(intersect( treePolygonTemp, forestPolyshape' num2str(v) '));']);
                        [~, forestTypeInd ] = max(forestCount);%it choose the first by default if more maxs
                    end%forest type
                    if sum(forestCount) == 0 %if points do not belonging to a known forest type
                        forestTypeInd = 1; %generic = populus
                    end

                    %estimate biomass of individual trees through species-dependent
                    %allometric relationships
                    treeHeight  = max(ptCloud(cellIndex,3)) - min(ptCloud(cellIndex,3)); %m
                    treeDiameter  = 0.01 * (treeHeight/aParameter(forestTypeInd))^(1/bParameter(forestTypeInd)); %m
                    treeBiomass  = (0.25 * pi * treeDiameter^2) * treeHeight * greenWoodDensity(forestTypeInd); %kg

                    %estimate the associated error
                    treeVariance = (treeBiomass^2)*( (treeDensityStd/greenWoodDensity(forestTypeInd))^2 + (-2*AStd(forestTypeInd)/(aParameter(forestTypeInd)*bParameter(forestTypeInd)))^2 + (-2*log(treeHeight/aParameter(forestTypeInd))*BStd(forestTypeInd)/(bParameter(forestTypeInd)^2))^2 + ((bParameter(forestTypeInd)+2)*HStd/(bParameter(forestTypeInd)*treeHeight))^2 );
                    treeStd = sqrt(treeVariance);

                    %write json file
                    listingtrees = cell({j});
                    filenameTrees = strcat('T_',string(cell2mat(listingtrees)));
                    jsonBiom = jsonUpdate(jsonBiom, 'tree', filenameTrees, treeBiomass,polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                    jsonErr = jsonUpdate(jsonErr, 'tree', filenameTrees, treeStd,polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                    if writeTotalJson
                        jsonBiomTot = jsonUpdate(jsonBiomTot, 'tree', filenameTrees, treeBiomass, polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                        jsonErrTot = jsonUpdate(jsonErrTot, 'tree', filenameTrees, treeStd, polyshape(ptCloudRoi(treeBoundaryVertices,1), ptCloudRoi(treeBoundaryVertices,2)));
                    end

                end

            end%loop over trees

            ptCloud(cellIndex,1:2) = NaN;

        end %if islands intersect the tile

        %remove point representing trees inside SA from the cloud
        nanIndex = find(isnan(ptCloud(:,1)));
        ptCloud(nanIndex,:) = [];

        %message
        fprintf('Tree biomass : computed\n')

        %%%%update vegetation mask
        vegetationMask(vegetationMask==-1) = 0;

        %message
        fprintf('%s\n', 'Vegetation mask : updated')


        %%%%detect individual shrubs
        if ~isempty(saIntersect.Vertices)
            fileOngoing = num2str(cell2mat(cell({'Area_',rgbFileName})));
            [~,fileOngoingName,~] = fileparts(num2str(cell2mat(cell({'Area_',rgbFileName}))));

            if ~existingClippedClouds

                %message
                fprintf('%s\n', 'Clipping shrub point clouds')

                %detect clusters within the vegetation mask
                clusters = bwboundaries(vegetationMask, 8,'noholes');

                %remove clusters that are too small to be vegetation
                [clusterRows,~] = cellfun(@size, clusters, 'UniformOutput', false);%it returns number of cells
                indToRemove = find( cell2mat(clusterRows) <= clusterThreshold);
                clusters(indToRemove) = [];
                xArray = tileX(1,:);
                yArray = tileY(:,1);
                clusterNumber = size(clusters, 1);

                %create polygons from clusters
                clusterPolygon = polyshape(clusterNumber,1);
                parfor k = 1 : clusterNumber
                    warning('off', 'all')
                    thisCluster = clusters{k};
                    xCluster = xArray(thisCluster(:,2));
                    yCluster = yArray(thisCluster(:,1));
                    clusterPolygon(k,1) = polyshape(xCluster',yCluster);
                end

                %separate multiple polygons
                shrubMask = polyshape();
                i = 0;
                count= 0;
                for k = 1 : size(clusterPolygon,1)
                    if ~isempty(clusterPolygon(k).Vertices)
                        inds = find(isnan(clusterPolygon(k).Vertices(:,1)));
                        if isempty(inds)
                            i = i+1;
                            shrubMask(i,:) = polyshape( clusterPolygon(k).Vertices(:,1)', clusterPolygon(k).Vertices(:,2)');
                        else
                            indStart = 1;
                            for j = 1 : numel(inds)+1
                                if j > numel(inds)
                                    indEnd = numel(clusterPolygon(k).Vertices(:,1));
                                    i = i+1;
                                    shrubMask(i,:) = polyshape( clusterPolygon(k).Vertices(indStart:indEnd,1)', clusterPolygon(k).Vertices(indStart:indEnd,2)');
                                else
                                    indEnd = inds(j)-1;
                                    i = i+1;
                                    shrubMask(i,:) = polyshape( clusterPolygon(k).Vertices(indStart:indEnd,1)', clusterPolygon(k).Vertices(indStart:indEnd,2)');
                                    indStart = inds(j)+1;
                                end
                            end
                        end
                    end
                end

                %optional very time consuming
                %         %compute the centroid and area of each polygon
                %         centroids = zeros(clusterNumber,4);
                %         [centroids(:,1),centroids(:,2)] = centroid(clusterPolygon(:,1));
                %         centroids(:,3) = area(clusterPolygon(:,1));
                %
                %         %nan removal
                %         nanInd = isnan(centroids);
                %         centroids((find(nanInd(:,1) == 1)),:) = [];
                %         clusterPolygon((find(nanInd(:,1) == 1)),:) = [];
                %
                %         %order centroids from the largest area to the smallest one
                %         centroids(:,4) = (1 : 1 : size(centroids,1)); %ids
                %         centroids = sortrows(centroids,3,'descend');
                %         clusterPolygon = clusterPolygon(centroids(:,4),:);
                %
                %         %buffer around the smallest clusters
                %         parfor k = 1 : size(shrubsPolygon,1) %par
                %             warning('off', 'all')
                %             if centroids(k,3) <= pi*(clusterRadius^2)
                %                 clusterPolygon(k,1) = polybuffer(centroids(k,1:2),'line',clusterRadius);
                %             end
                %         end
                %
                %         %define mask corresponding to the detected shrubs
                %         shrubMask = union(shrubsPolygon);


                %message
                fprintf('%s\n', 'Shrub masks: defined')


                %%%%clip the point cloud according to the detected shrubs

                %clip the point cloud
                shrubNumber = size(shrubMask,1);
                ptCloudX = ptCloud(:,1);
                ptCloudY = ptCloud(:,2);
                ptCloudZ = ptCloud(:,3);
                parfor t = 1 : shrubNumber
                    warning('off', 'all')
                    [maskX, maskY] = poly2cw(shrubMask(t,1).Vertices(:,1), shrubMask(t,1).Vertices(:,2));
                    clusterPolygon = [maskX maskY];
                    pointCloudInIndex = inpolygon(ptCloudX, ptCloudY, maskX, maskY);
                    Xp = ptCloudX(pointCloudInIndex == 1,:);
                    Yp = ptCloudY(pointCloudInIndex == 1,:);
                    Zp = ptCloudZ(pointCloudInIndex == 1,:);

                    if size(Xp,1) >= 4 %otherwise a volume cannot be defined
                        listingShrubs = cell({t});
                        filenameShrubs = cell2mat(listingShrubs);

                        %create local output directory
                        if ~exist(strcat(outputFolder, '/Shrubs/', fileOngoingName, '/'), 'dir')
                            mkdir(strcat(outputFolder, '/Shrubs/', fileOngoingName, '/'))
                        end

                        %save clipped cloud
                        writetable(table(Xp,Yp,Zp),strcat(outputFolder, '/Shrubs/',fileOngoingName,'/S_',num2str(filenameShrubs),'.txt'))

                    end

                end %shrubNumber

                %message
                fprintf('%s\n', 'Shrub point clouds: clipped')

            end%if clipped clouds do not exist

            %%%compute shrub volume and biomass

            %message
            fprintf('%s\n', 'Computing shrub biomass')

            %retrieve shrub list
            shrubList = strcat(outputFolder,'/Shrubs/',fileOngoingName,'/*.txt');
            shrubNames = struct2cell(dir(shrubList));
            namesArray = sort(shrubNames(1,:));
            shrubNumber = size(shrubNames,2);

            %initialization
            shrubVolume = NaN(shrubNumber,1);
            shrubZMax = NaN (shrubNumber,1);
            shrubZMin = NaN (shrubNumber,1);
            shrubNPoints = NaN (shrubNumber,1);
            shrubNormalizedPoints = NaN (shrubNumber,1);
            shrubAreaMax = NaN (shrubNumber,1);
            shrubBiomass = NaN (shrubNumber,1);
            shrubVariance = NaN(shrubNumber,1);
            shrubStd = NaN(shrubNumber,1);
            maxAreaX = []; maxAreaY = [];

            %loop over shrubs
            for t = 1 : shrubNumber

                fprintf('shrub #%s of %s\n', num2str(t), num2str(shrubNumber))

                %import
                inputShrub = readtable((strcat(outputFolder, '/Shrubs/',fileOngoingName,'/',cell2mat(namesArray(1,t)))));
                inputArray = table2array(inputShrub);

                %remove equal points
                [~,idx]= unique(inputArray ,'rows','first');
                inputArray  = inputArray(sort(idx),:);
                z = inputArray(:,3);
                zMin = min(z);
                zMax = max(z);
                maxAreaX = []; maxAreaY = [];

                %remove ground points
                if size(find((z-zMin) > groundThreshold),1) > 5
                    inputArrayIndex = find(z > zMin + groundThreshold);
                    inputArray = inputArray(inputArrayIndex,:);
                    x = inputArray(:,1);
                    y = inputArray(:,2);
                    z = inputArray(:,3);
                    shrubNPoints(t,1) = size(z,1);

                    if size(z,1) > 0  %not only ground

                        %height classes
                        shrubZMin(t,1) = min(z);
                        shrubZMax(t,1) = max(z);
                        zClasses = shrubZMin(t,1) : dh : shrubZMax(t,1)+dh;
                        classesNumber = size(zClasses,2)-1;

                        if shrubNPoints(t,1) < 5 %convex hull slicing method

                            cloud = [inputArray(:,1) inputArray(:,2) inputArray(:,3)];
                            [~, shrubVolume(t,1)] = boundary(cloud);
                            [k, shrubAreaMax(t,1)] = boundary(cloud(:,1:2),0);
                            maxAreaX = cloud(k,1);
                            maxAreaY = cloud(k,1);
                            shrubAreaMax(t,1) = polyarea(maxAreaX, maxAreaY);
                            shrubVolume(t,1) = shrubVolume(t,1)*volumeCorrection;%m3 corrected!!!

                        else %concave hull slicing method

                            %dh array
                            dhArray0 = ones(1,classesNumber) * dh;

                            %initial classes
                            zClasses0 = discretize(z,zClasses);
                            zClasses0Size  = size(zClasses,2)-1;
                            pointPerClass0 = zeros(1,zClasses0Size);
                            for i = 1:zClasses0Size
                                pointPerClass0(1,i) = length(find(zClasses0==i));
                            end

                            %merge classes with too few points starting from the ground
                            zClasses1 = zClasses0;
                            dhArray1 = dhArray0;
                            i = 1;
                            j = 1;
                            while i <=  zClasses0Size
                                pointsInIndex = find(zClasses1==i);
                                pointsInNumber = size(pointsInIndex);
                                zClasses1(pointsInIndex,:) = (j).*ones(pointsInNumber);
                                dhArray1(1,j) = dhArray0(1,i);
                                %if the upper layer
                                if i  ==  zClasses0Size
                                    %merge the upper layer with the one below
                                    if pointsInNumber(1) > 1 && pointsInNumber(1) < pointThreshold
                                        zClasses1(pointsInIndex,:) = (j-1).*ones(pointsInNumber);
                                        dhArray1(1,j-1) = dhArray1(1,j)+dhArray1(1,j-1);
                                        break
                                    else
                                        zClasses1(pointsInIndex,:) = (j).*ones(pointsInNumber);
                                        dhArray1(1,j) = dhArray0(1,i);
                                        j=j+1;
                                        break
                                    end
                                    %if not the upper layer
                                elseif i  <  zClasses0Size
                                    %too few points
                                    if pointsInNumber(1) < pointThreshold
                                        pointsInNumberJ = [0,0];
                                        %merge layers until the number of points per layer is
                                        %adequate
                                        while pointsInNumberJ(1) < pointThreshold && i <  zClasses0Size
                                            pointsInIndexNew = find(zClasses1==i+1);
                                            zClasses1(pointsInIndexNew,:) = (j).*ones(size(pointsInIndexNew));
                                            pointsInIndexJ = find(zClasses1==j);
                                            pointsInNumberJ = size(pointsInIndexJ);
                                            dhArray1(1,j) = dhArray1(1,j)+dhArray0(1,i+1);
                                            %dhArray1(1,j+1) = dhArray1(1,j+1)-dhArray0(1,i+1);
                                            if pointsInNumberJ < pointThreshold
                                                i = i + 1;
                                            else
                                                i = i + 2;
                                            end
                                        end%while
                                        %enough points
                                    else
                                        zClasses1(pointsInIndex,:) = (j).*ones(pointsInNumber);
                                        dhArray1(1,j) = dhArray0(1,i);
                                        i = i+1;
                                    end
                                    j = j+1;
                                    %if over the upper layer
                                else
                                    break
                                end
                            end%while

                            %update the class size
                            dhArray1 = dhArray1(:,1:j-1);
                            zClasses1Size = size(dhArray1,2);
                            pointPerClass1 = zeros(1,zClasses1Size);
                            for i = 1:zClasses1Size
                                pointPerClass1(1,i) = length(find(zClasses1==i));
                            end
                            if sum(pointPerClass0) ~= sum(pointPerClass1)
                                error('wrong merging')
                            end

                            %categorize the areas of the layers
                            areas = zeros(1,zClasses1Size);
                            for i =  1 :zClasses1Size
                                hull = EBE_concave_hull(table(x(zClasses1 == i),y(zClasses1 == i)),7);
                                areas(i) = polyarea(hull.x, hull.y);
                            end
                            layerAreaMean = mean(areas);
                            layerAreaStd = std(areas);
                            areaCategories = zeros(1,zClasses1Size);
                            for i = 1 : zClasses1Size
                                if areas(i) - layerAreaMean >=  0
                                    areaCategories(i) = fix((areas(i)-layerAreaMean)/layerAreaStd)+1;
                                else
                                    areaCategories(i) = fix((areas(i)-layerAreaMean)/layerAreaStd)-1;
                                end
                            end

                            %merge adjacent layers with similar areas
                            zClasses2 = zClasses1;
                            dhArray2 = dhArray1.*0;
                            i = 1;
                            j = 1;
                            while i <= zClasses1Size
                                if i == zClasses1Size
                                    pointsInIndexNew = find(zClasses2==i);
                                    zClasses2(pointsInIndexNew,:) = (j).*ones(size(pointsInIndexNew));
                                    dhArray2(1,j) = dhArray1(1,i);
                                    j = j+1;
                                    break
                                end
                                pointsInIndexNew = find(zClasses2==i);
                                zClasses2(pointsInIndexNew,:) = (j).*ones(size(pointsInIndexNew));
                                dhArray2(1,j) = dhArray1(1,i);
                                if areaCategories(i)  ==  areaCategories(i+1)
                                    %merge
                                    while areaCategories(i)  ==  areaCategories(i+1)
                                        pointsInIndexNew = find(zClasses2==i+1);
                                        zClasses2(pointsInIndexNew,:) = (j).*ones(size(pointsInIndexNew));
                                        dhArray2(1,j) = dhArray2(1,j) + dhArray1(1,i+1);
                                        i = i+1;
                                        if i  ==  zClasses1Size
                                            break
                                        end
                                    end
                                    i = i + 1;

                                else
                                    pointsInIndexNew = find(zClasses2==i);
                                    zClasses2(pointsInIndexNew,:) = (j).*ones(size(pointsInIndexNew));
                                    dhArray2(1,j) = dhArray1(1,i);
                                    i = i+1;
                                end
                                j = j+1;
                            end

                            %update class size
                            dhArray2 = dhArray2(:,1:j-1);
                            zClasses2Size = size(dhArray2,2);
                            pointPerClass2 = zeros(1,zClasses2Size);
                            for i = 1:zClasses2Size
                                pointPerClass2(1,i) = length(find(zClasses2==i));
                            end
                            if sum(pointPerClass2) ~= sum(pointPerClass1)
                                error('wrong merging of layers with similar areas')
                            end

                            %compute areas and volumes
                            areas = zeros(1,zClasses2Size);
                            for i =  1 :zClasses2Size
                                hull = EBE_concave_hull(table(x(zClasses2 == i),y(zClasses2 == i)),7);
                                areas(i) = polyarea(hull.x, hull.y);
                                if areas(i) > polyarea(maxAreaX, maxAreaY)
                                    maxAreaX = hull.x;
                                    maxAreaY = hull.y;
                                    shrubAreaMax(t,1) = polyarea(maxAreaX, maxAreaY);
                                end
                            end
                            layerAreaMean = mean(areas);
                            layerAreaStd = std(areas);
                            layerVolume = areas.*dhArray2;
                            shrubVolume(t,1) = sum(layerVolume)*volumeCorrection;%m3 corrected!!!

                        end

                        if  ~isnan(shrubVolume(t,1)) && shrubVolume(t,1) > 0 && ~isempty(maxAreaX) && ~isempty(maxAreaY)

                            %remove volume due to trees enclosed within the
                            %shrub cluster
                            shrubPolygon = polyshape(maxAreaX,maxAreaY);
                            for i = 1:size(treeOver ,1)
                                treeOverlapPolygon = intersect(treePolygon(i,:),shrubPolygon);
                                overlapArea = area(treeOverlapPolygon);
                                if overlapArea >0
                                    shrubVolume(t,1) = shrubVolume(t,1) * (1- overlapArea/area(shrubPolygon));
                                end
                            end

                            %biomass computation
                            shrubBiomass(t,1) = shrubVolume(t,1) * bulkDensity;%kg

                            %compute the associated error
                            shrubVariance(t,1) = (shrubDensityStd^2)*(ALSVolStd^2) + (shrubDensityStd^2)*(shrubVolume(t,1)^2) + (ALSVolStd^2)*(bulkDensity^2);
                            shrubStd(t,1) = sqrt(shrubVariance(t,1));

                            %write json file
                            [~,filenameShrubs,~] = fileparts(string(cell2mat(namesArray(1,t))));
                            nanvalue = find(isnan(shrubPolygon.Vertices(:,1)));
                            if isempty(nanvalue)
                                jsonBiom = jsonUpdate(jsonBiom, 'shrub', filenameShrubs, shrubBiomass(t,1),shrubPolygon);
                                jsonErr = jsonUpdate(jsonErr, 'shrub', filenameShrubs, shrubStd(t,1),shrubPolygon);
                                if writeTotalJson
                                    jsonBiomTot = jsonUpdate(jsonBiomTot, 'shrub', filenameShrubs, shrubBiomass(t,1),shrubPolygon);
                                    jsonErrTot = jsonUpdate(jsonErrTot, 'shrub', filenameShrubs, shrubStd(t,1),shrubPolygon);
                                end
                            else %multypolygon if it contains at least 2 polygons
                                shrubPolygonCorr = shrubPolygon;
                                shrubPolygonCorr.Vertices = [];
                                shrubPolygonCorr.Vertices = [shrubPolygonCorr.Vertices; shrubPolygon.Vertices(1:nanvalue(1)-1,:)];
                                for i = 1:size(nanvalue,1)
                                    indstart = nanvalue(i) +1;
                                    if size(nanvalue,1) == 1
                                        indEnd = size(shrubPolygon.Vertices,1);
                                    end
                                    if size(nanvalue,1) ~= 1 && i < size(nanvalue,1)
                                        indEnd = nanvalue(i+1)-1;
                                    elseif size(nanvalue,1) ~= 1 && i == size(nanvalue,1)
                                        indEnd = size(shrubPolygon.Vertices,1);
                                    end
                                    jsonBiom = jsonUpdate(jsonBiom, 'shrub', filenameShrubs, shrubBiomass(t,1),shrubPolygonCorr);
                                    jsonErr = jsonUpdate(jsonErr,'shrub', filenameShrubs, shrubStd(t,1),shrubPolygonCorr);
                                    if writeTotalJson
                                       jsonBiomTot = jsonUpdate(jsonBiomTot, 'shrub', filenameShrubs, shrubBiomass(t,1),shrubPolygonCorr);
                                       jsonErrToT = jsonUpdate(jsonErrTot,'shrub', filenameShrubs, shrubStd(t,1),shrubPolygonCorr); 
                                    end
                                end
                            end


                        end

                    end %if not only ground

                end %remove ground

            end %shrubNumber

        end%if the SA intersects the  tile
      
      
        %%%%save json files
        encodedBiom = jsonencode(jsonBiom);
        fid = fopen(strcat(outputFolder,'/exportedJson/Single/', fileName,'_biomass','.geojson'),'w');
        fprintf(fid,encodedBiom);
        fclose(fid);
        encodedErr = jsonencode(jsonErr);
        fid = fopen(strcat(outputFolder,'/exportedJson/Single/', fileName,'_error','.geojson'),'w');
        fprintf(fid,encodedErr);
        fclose(fid);

        %message
        fprintf('%s\n', 'Shrub biomass : computed')

        %delete variables
        shrubMask = polyshape();
        clear ptCloud

    end % if RoI overlaps the tile

end %loop over input areas


%%%%%%%%%%%%%%% output generation %%%%%%%%%%%%%%%

%message
fprintf('%s\n', 'Generating output files')

%%%% save json files
if writeTotalJson
    encodedBiomTot = jsonencode(jsonBiomTot);
    fid = fopen(strcat(outputFolder,'/exportedJson/', jsonFileName, '_biomass','.geojson'),'w');
    fprintf(fid,encodedBiomTot);
    fclose(fid);
    encodedErrTot = jsonencode(jsonErrTot);
    fid = fopen(strcat(outputFolder,'/exportedJson/', jsonFileName, '_error','.geojson'),'w');
    fprintf(fid,encodedErrTot);
    fclose(fid);
end

%message
fprintf('%s\n', 'END')
