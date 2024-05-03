%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Biomass Estimation EBE - version 1.0
%
% Individual tree detection and measurement (ITDM) module
%
% Melissa Latella
%
% Created: September 2020
% Last update: February 2021
%
% melissa.latella@cmcc.it
%
% Cite as: "Latella, M., Sola, F., & Camporeale, C. (2021). A Density-Based
%          Algorithm for the Detection of Individual Trees from LiDAR Data.
%           Remote Sensing, 13(2), 322
%           and "Latella, M., Raimondo, T., Belcore, L., Salerno, L.,
%                 and Camporeale, C. (2022), On the integration of LiDAR
%                   and field data for riparian biomass estimation".
%
% Aim: To identify and measure individual trees within a stand by
%       analysing 3D point clouds.
%
% Method: the input point cloud must have the z coordinates expressed as a
% relative height (i.e. absolute level - ground level).
% It is preferable to use "all return" normalised LiDAR point clouds as
% input.
% The user can choose to identify the stems or the tree tops.
% -> Stem identification: the algorithm first removes the outliers
% according to a user-defined height threshold. It also removes all the
% points below 1.4 m, which represent shrubs, brushes and grasses. Then,
% for each point, it clips the cloud and computes the most frequent radius
% within which the other points of the clouds can be found. It computes the
% point density in the area of radius equal the most frequent radius and
% identifies the stem position as the location of the maximum point density
% within the area defined by the most frequent radius.
% It finally filters the double-counted stems.
% -> Tree top identification: same as above but the algorithm then shifts
% the coordinates at the location of the closest height maximum around the
% stems so to identify the top.
% Optionally, it filters the anthropic objects in the stands after that
% the user has defined their range of height.
%
% Dependencies:
% 1. Parallel Toolbox  by Matlab (otherwise substitute "parfor" with "for"
%    in the script and comment rows 91-92).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%%%%%%%%%%%%%%% parameter setting %%%%%%%%%%%%%%%

%number of cores for parallelization
coreNumber = 8;

%x input
xColumn = 1;

%y input
yColumn = 2;

%height input
heightColumn = 3;

%clipping radius
clippingRadius = 20; %m

%minimum stem spacing
%used to detect double-counted trees (9999 means automated calculation)
stemSpacingPreference = 9999;

%above threshold outliers (birds, errors)
aboveOutlierThreshold = 50; %m

%below threshold outliers (bushes, shrubs)
belowOutlierThreshold = 1.4; %m

%thresholds for local outliers due to anthropic objects
anthropicRangeNumber = 0; %number of ranges
anthropicLowerThreshold = NaN(anthropicRangeNumber,1);%array with the lowest values of the ranges
anthropicUpperThreshold = NaN(anthropicRangeNumber,1);%array with the upper values of the ranges

%output preference: false = only stems, true = stems + treetops
treetopPreference = true;

%%%%%%%%%%%%%%% parallel definition %%%%%%%%%%%%%%%

delete(gcp('nocreate')) %it requires the Parallel computing toolbox
poolobj = parpool('local',coreNumber);

%%%%%%%%%%%%%%% selection of the point cloud files %%%%%%%%%%%%%%%

%%%%folder definition
inputFolder = 'xxx';
inputFormat = '.las';%.csv or .txt are allowed
inputAOI = shaperead('xxx.shp');
[AOIx,AOIy] = poly2cw(inputAOI.X(:),inputAOI.Y(:));

%%%%file list creation
listing = struct2cell(dir(strcat(inputFolder,'\*',inputFormat)));
namesArrayTxt = sort(listing(1,:));
%number of files
T = size(listing,2);

%%%%output directory
currentFolder = pwd;
if ~exist('Output', 'dir')
    mkdir('Output')
end


%%%%%%%%%%%%%%% loop over the file list %%%%%%%%%%%%%%%

for t = 1 : T

    %reading input file
    pointCloudFile = cell2mat(namesArrayTxt(1,t))
    lasReader = lasFileReader(strcat(inputFolder,pointCloudFile));
    [ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","ScanAngle");
    pointCloud = ptCloud.Location;
    pointCloudFile(find(pointCloudFile == '.', 1 ,'last') : end) = [];


    %%%%%%%%%%%%%%% cloud pre-processing %%%%%%%%%%%%%%%

    %%%%removal of outliers
    %above threshold outliers (birds, errors)
    pointCloud(pointCloud(:, heightColumn) > aboveOutlierThreshold , :) = [];
    %below threshold (shrubs, bushes, grasses)
    pointCloud(pointCloud(:, heightColumn) < belowOutlierThreshold , :) = [];

    if ~isempty(pointCloud)
        %%%%%%%%%%%%%%% stem detection %%%%%%%%%%%%%%%

        %%%%point coordinates
        x = pointCloud(:,xColumn);
        y = pointCloud(:,yColumn);
        z = pointCloud(:,heightColumn);

        %%%%AoI
        ind = find(inpolygon(x, y, AOIx, AOIy));

        if ~isempty(ind)
              
            x = x(ind);
            y = y(ind);
            z = z(ind);     

            %%%%initial point number
            pointNumber = length(x);

            %%%%variable initialization
            mostFrequentRadius = rand(pointNumber,1);
            density = NaN(pointNumber,1);
            stemX = NaN(pointNumber,1);
            stemY = NaN(pointNumber,1);
            stemHeight = NaN(pointNumber,1);
            stemMostFrequentRadius = NaN(pointNumber,1);

            %%%%most frequent radius computation and point density around each point
            parfor i = 1 : pointNumber
                distance  =  sqrt((x-x(i,1)).^2 + (y-y(i,1)).^2);
                rangeIndex  =  find(distance < clippingRadius);
                if length(rangeIndex) < 2%if not enough point are found within the clipping radius
                    mostFrequentRadius(i,1) = clippingRadius;
                else
                    mostFrequentRadius(i,1) = mostFreqRadiusComputation(distance(rangeIndex));
                end
                rangeIndexUpdated  =  find(x(i,1) > (x - mostFrequentRadius(i,1)) & x(i,1) < (x + mostFrequentRadius(i,1)) & y(i,1) > (y - mostFrequentRadius(i,1)) & y(i,1) < (y + mostFrequentRadius(i,1)) );
                density(i,1)  =  length(rangeIndexUpdated)/(4*mostFrequentRadius(i,1)^2);%square of side equal to the mostFrequentRadius
            end

            %%%%stem position
            parfor i = 1 : pointNumber
                rangeIndex  =  find(x(i,1) > (x - mostFrequentRadius(i,1)) & x(i,1) < (x + mostFrequentRadius(i,1)) & y(i,1) > (y - mostFrequentRadius(i,1)) & y(i,1) < (y + mostFrequentRadius(i,1)));
                rangeDensity = zeros(pointNumber,1); %density of the points inside the area defined by the most frequent radius
                rangeDensity(rangeIndex)  =  density(rangeIndex);

                if rangeDensity(i) == max(rangeDensity) % if the i-th point has the maximum density in its around, it is a stem
                    stemX(i,1) = x(i,1);
                    stemY(i,1) = y(i,1);
                    stemHeight(i,1) = z(i,1);
                    stemMostFrequentRadius(i,1) = mostFrequentRadius(i,1);
                end

            end
            nanIndex = find(isnan(stemX));
            stemX(nanIndex) = []; stemY(nanIndex) = [];
            stemHeight(nanIndex) = [];
            stemMostFrequentRadius(nanIndex) = [];

            %%%%saving arrays in new variables
            stemX0 = stemX;
            stemY0 = stemY;
            stemHeight0 = stemHeight;
            stemMostFrequentRadius0 = stemMostFrequentRadius;

            %%%%output
            treeListS0 = table(stemX0, stemY0); %coordinates of stems
            writetable(treeListS0, strcat(currentFolder,'\Output\',pointCloudFile,'_stems_beforeFilter.txt'));

            %%%%most frequent stem spacing

            stemX = stemX0;
            stemY = stemY0;
            stemHeight = stemHeight0;
            stemMostFrequentRadius = stemMostFrequentRadius0;
            stemNumber = length(stemX);
            if stemSpacingPreference == 9999 %preference set at the beginning
                stemSpacing = rand(1,1); %initialization
                stemSpacing = stemSpacingComputation(stemSpacing, stemX,stemY,stemNumber);
            else %if the user has specified a site-specific stem spacing
                stemSpacing = stemSpacingPreference;
            end

            %%%%Filter to remove double-counted stems
            for i  =  1 : stemNumber
                if stemX(i) ~= 9999
                    doubleIndex =  find((((stemX - stemX(i)).^2 + (stemY - stemY(i)).^2).^0.5) <= stemSpacing);
                    if length(doubleIndex) > 1 %it counts also the point itself
                        [~, newIndex]  = max(stemHeight(doubleIndex)); %if there are more points of the same height it arbitrary chooses the first one
                        %newPosition = doubleIndex(newIndex);
                        %removal of the chosen "true" tree top from the double-counted list
                        doubleIndex(newIndex) = [];
                        %label for the "false" tree top not to consider in the following iterations
                        if ~isempty(doubleIndex)
                            stemX(doubleIndex) = 9999;
                        end
                    end
                end
            end
            %removal of all false-labelled stems
            falseIndex = find(stemX == 9999);
            stemX(falseIndex) = [];
            stemY(falseIndex) = [];
            stemMostFrequentRadius(falseIndex) = [];
            falseNumber = length(falseIndex);
            stemNumber  =  length(stemX);
            treeNumber = stemNumber;

            %%%%write results on table
            treeListS = table(stemX, stemY); %coordinates of stems

            %%%%write output file
            writetable(treeListS, strcat(currentFolder,'\Output\',pointCloudFile,'_stems_SS', string(round(stemSpacing.*10,0)) ,'.txt'));


            %%%%%%%%%%%%%%% tree top detection %%%%%%%%%%%%%%%

            if treetopPreference

                %%%%variable initialization
                maxHeight = NaN(stemNumber,1);
                topX = NaN(stemNumber,1);
                topY = NaN(stemNumber,1);
                topHeight = NaN(stemNumber,1);

                %%%%tree top position
                [topX, topY, topHeight] = findTopPosition(topX, topY, topHeight, stemX,stemY,stemMostFrequentRadius, stemNumber,x,y,z,pointNumber);

                %%%%write results on table
                treeListT = table(topX, topY, topHeight);% x,y coordinates, and height of treetops

                %%%%write output file
                writetable(treeListT, strcat(currentFolder,'\Output\',pointCloudFile,'_tops_SS',  string(round(stemSpacing.*10,0)) ,'_beforeAnthropicFilter.txt'));

                %%%%%%%%%%%%%%% removal of anthropic objects %%%%%%%%%%%%%%%

                if anthropicRangeNumber > 0 %preference set at the beginning

                    aroundRadius = 20;
                    localOutlierNumber = 0;

                    for i = 1 : anthropicRangeNumber

                        anthropicIndex = find(topHeight >= anthropicLowerThreshold(i) & topHeight <= anthropicUpperThreshold(i));
                        if ~isempty(anthropicIndex)
                            anthropicNumber = length(anthropicIndex);
                            for j = 1 : anthropicNumber
                                rangeIndex = find(topX(anthropicIndex(j)) > (topX - aroundRadius) & topX(anthropicIndex(j)) < (topX + aroundRadius) & topY(anthropicIndex(j)) > (topY - aroundRadius) & topY(anthropicIndex(j)) < (topY + aroundRadius));
                                %exclude labelled with 9999
                                labelledIndex = find(topX(anthropicIndex(j) == 9999));
                                rangeIndex(labelledIndex) = [];
                                if topHeight(anthropicIndex(j)) < (mean(topHeight(rangeIndex), 'omitnan') - 2*std(topHeight(rangeIndex), 'omitnan')) || topHeight(anthropicIndex(j)) > (mean(topHeight(rangeIndex), 'omitnan') + 2*std(topHeight(rangeIndex), 'omitnan'))
                                    %label for the "outlier" tree top not to consider in the following iterations
                                    topX(anthropicIndex(j))
                                end
                            end
                        end
                        %removal of all outlier-labelled tree tops
                        localOutlierIndex = find(topX == 9999);
                        topX(localOutlierIndex) = [];
                        topY(localOutlierIndex) = [];
                        topHeight(localOutlierIndex) = [];
                        localOutlierNumber = localOutlierNumber + length(localOutlierIndex);

                    end%i

                    %%%%write results on table
                    treeListT = table(topX, topY, topHeight);% x,y coordinates, and height of treetops

                    %%%%write output file
                    writetable(treeListT, strcat(currentFolder,'\Output\',pointCloudFile,'_tops_SS',  string(round(stemSpacing.*10,0)) ,'.txt'));

                end %if anthropic objects

                treeNumber = length(topHeight);

            end %if "tree top detection" preference

            %%%%store information file
            stemSpacingValue = stemSpacing;
            countedTrees = treeNumber;

            %%%%write information file
            summary = table(stemSpacingValue, countedTrees);
            writetable(summary, strcat(currentFolder,'\Output\',pointCloudFile,'_StemSpacing_results.txt'))

        end%if isempty
    end%if isempty
end%t

%%%%Closing parallel pool
delete(poolobj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%function to compute the most frequent radius
function mostFrequentRadius = mostFreqRadiusComputation(distance)

%1/distance array
inverseDistance = 1./distance;
inverseDistanceOrdered = sort(inverseDistance);
inverseDistanceOrdered(isinf(inverseDistanceOrdered)) = []; %to eliminate Inf value

%options
angularFrequency = 2*pi;
nfft = size(inverseDistanceOrdered,1);
kaiserWindow = ones(nfft,1);

%power spectrum
distanceFFT = fft(inverseDistanceOrdered,nfft);
angularFrequencyStep = angularFrequency /nfft;
frequency1 = angularFrequencyStep*(0:nfft-1);

%fix error
NyquistFrequency = angularFrequency/2;
angularFrequencyResolutionHalf = angularFrequencyStep/2;
if rem(nfft,2) == 1 %if odd
    halfNfft = (nfft+1)/2;
    % Adjust points on either side of NyquistFrequency
    frequency1(halfNfft)   = NyquistFrequency - angularFrequencyResolutionHalf;
    frequency1(halfNfft+1) = NyquistFrequency + angularFrequencyResolutionHalf;
else %if even
    halfNfft = (nfft/2)+1;
    frequency1(halfNfft) = NyquistFrequency;
end

%auto spectrum
frequency1(nfft) = angularFrequency - angularFrequencyStep;
windowCompensation = kaiserWindow'*kaiserWindow;
autoSpectrum1 = distanceFFT.*conj(distanceFFT)/windowCompensation;
if rem(nfft,2) == 1 %if odd
    rangeSelection = 1:(nfft+1)/2;
    autoSpectrumUnscaled = autoSpectrum1(rangeSelection,:); % to take just [0,pi]
    autoSpectrum = [autoSpectrumUnscaled(1,:); 2*autoSpectrumUnscaled(2:end,:)];
else %if even
    rangeSelection = 1:nfft/2+1;
    autoSpectrumUnscaled = autoSpectrum1(rangeSelection,:); % to take just [0,pi]
    autoSpectrum = [autoSpectrumUnscaled(1,:); 2*autoSpectrumUnscaled(2:end-1,:); autoSpectrumUnscaled(end,:)];
end
frequency = frequency1(rangeSelection);

%one-sided PSD [Power/freq]
powerSpectralDensity = autoSpectrum./angularFrequency;

%width of each frequency bin
if numel(frequency)>2% ~isempty(frequency)
    binWidth = NaN(length(frequency),1);
    binWidth(1:end-1,1) = diff(frequency);
    binWidth(end) = binWidth(end-1);

    %power within each bin and total power
    binPower = binWidth .* powerSpectralDensity;
    totalPower = sum(binPower(:));

    %central moment within the range
    meanFrequency = sum(binPower(:) .* frequency(:)) ./ totalPower;

    %most frequent radius
    mostFrequentRadius = 1 / meanFrequency;


else
    %most frequent radius
    mostFrequentRadius = 0;
end


end

%%%%function to compute the stem spacing according to the cloud features
function stemSpacing = stemSpacingComputation(stemSpacing,stemX,stemY,stemNumber)
clippingRadius2 = 100;%by default
mostFrequentSpacing = NaN(stemNumber,1);
parfor i = 1 : stemNumber
    distance  =  sqrt((stemX-stemX(i,1)).^2 + (stemY-stemY(i,1)).^2);
    rangeIndex  =  find(distance < clippingRadius2);
    if length(rangeIndex) < 2
        mostFrequentSpacing(i,1) = clippingRadius2;
    else
        mostFrequentSpacing(i,1) = mostFreqRadiusComputation(distance(rangeIndex));
    end
end
%median stem spacing
stemSpacing = median(mostFrequentSpacing, 'omitnan');%default value chosen after a trial-and-error process
end


%%%%function to find the closest treetops with respect to the identified
%%%%stem
function  [topX, topY, topHeight] = findTopPosition(topX, topY, topHeight,stemX,stemY,stemMostFrequentRadius, stemNumber,x,y,z,pointNumber)
parfor i = 1 : stemNumber
    rangeIndex  =  find(stemX(i,1) > (x - stemMostFrequentRadius(i,1)) & stemX(i,1) < (x +stemMostFrequentRadius(i,1)) & stemY(i,1) > (y-stemMostFrequentRadius(i,1)) & stemY(i,1)< (y +stemMostFrequentRadius(i,1)));
    rangeHeight = zeros(pointNumber,1);
    rangeHeight(rangeIndex)  =  z(rangeIndex);
    maxHeight = max(rangeHeight);
    maxHeightPosition = find(rangeHeight == maxHeight); %if more than one, it moves the tree at the first location (arbitrary)
    topHeight(i,1) = z(maxHeightPosition(1));
    topX(i,1) = x(maxHeightPosition(1));
    topY(i,1) = y(maxHeightPosition(1));
end
end
