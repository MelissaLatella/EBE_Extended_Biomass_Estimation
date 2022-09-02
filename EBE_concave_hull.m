%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extended Biomass Estimation EBE - version 2.0
%
% Concave hull module
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
% Aim: To find optimal concave hull for a given set of points.
%
% Method: based on "Moreira, A., & Santos, M. Y. (2007). Concave hull:
%                       A k-nearest neighbours approach for the computation
%                       of the region occupied by a set of points,
%                       GRAPP 2007 - International Conference on Computer
%                       Graphics Theory and Application".
%         and "Yan, Z., Liu, R., Cheng, L., Zhou, X., Ruan, X., & Xiao, Y.
%               (2019). A concave hull methodology for calculating the
%               crown volume of individual trees based on vehicle-borne
%               LiDAR data. Remote Sensing, 11(6), 623."
%
% Warning: it requires the function selfintersect.m by Antoni J. Canós 
% (2022). Fast and Robust Self-Intersections 
% (https://www.mathworks.com/matlabcentral/fileexchange/13351-fast-and-robust-self-intersections), 
% MATLAB Central File Exchange. Retrieved July 29, 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function computedPolygon = EBE_concave_hull(points0,kk0)

%%%%%%%%%%%%%%% preliminary operations %%%%%%%%%%%%%%%


%%%%prepare point cloud

%rename columns
points0.Properties.VariableNames = [{'x'},{'y'}];

%remove equal points
[~,idx]= unique(points0,'rows','first');
points0 = points0(sort(idx),:);


%%%%other setting
kk0 = max(kk0,3);%initial neighborhood size


%%%%%%%%%%%%%%% concave hull %%%%%%%%%%%%%%%

if size(points0,1) <= 5 

    hullID = boundary(table2array(points0),0);
    hull = points0(hullID,:);

% elseif size(points0,1) >= 200
% 
%     hullID = boundary(table2array(points0),1);
%     hull = points0(hullID,:);

else

    kk0 = min(kk0,(size(points0,1)-1));
    allInside = false;

    %initialization
    NPoints = size(points0,1);
    distances = array2table(NaN(NPoints,1));
    angles = array2table(NaN(NPoints,1));
    points0 = [points0,distances];
    points0.Properties.VariableNames(3) = {'distances'};
    points0 = [points0,angles];
    points0.Properties.VariableNames(4) = {'angles'};
    hull = array2table(NaN(NPoints+1,4));%max size
    hull.Properties.VariableNames = {'x','y','distances','angles'};
    points = points0;

    %find the point with the lowest y coordinate to be the first point
    if size(find(points.y == min(points.y)),1) == 2

        ind = find(points.y == min(points.y));
        indMinY = ind(points.x(ind) == max(points.x(ind)));
        firstPoint = points(indMinY,:);

    elseif size(find(points.y == min(points.y)),1) > 3

        ind = find(points.y == min(points.y));
        indMinY = ind(points.x(ind) == max(points.x(ind)));
        firstPoint = points(indMinY,:);

    else

        [~,indMinY] = min(points.y);
        firstPoint = points((indMinY),:);
    end

    while ~allInside && kk0 <= size(points0,1)-1

        step = 2;
        hull = array2table(NaN(NPoints+1,4));%max size
        hull.Properties.VariableNames = {'x','y','distances','angles'};
        hull(1,:) = firstPoint;
        indHull = 2;
        points = points0;

        %move to current point and remove first point from dataset
        currentPoint = firstPoint;
        points(indMinY,:) = [];
        previousPoint = firstPoint;
        previousPoint.x = previousPoint.x-100;

        %increase neighborhood size
        kk0 = kk0+1;

        while ( ((currentPoint.x ~= firstPoint.x) && (currentPoint.y ~= firstPoint.y)) || step == 2) && size(points,1) > 0

            %                 if NPoints > 10 && step >= 8
            %                     points = [points;firstPoint];
            %                 elseif step >= 6
            %                     points = [points;firstPoint];
            %                 end
            if size(hull,1)>2
                points = [points;firstPoint];
            end

            if kk0 <= size(points0,1)
                kk = kk0;
                its = true;

                %find potential candidate points
                kkInd1 = kk;
                kkInd2 = size(points,1);
                if kkInd2<kk
                    kkInd1 = kkInd2;
                end

                for kk = kkInd1  : kkInd2

                    %select kk nearest points
                    distances = sqrt((currentPoint.y - points.y).^2 + (currentPoint.x - points.x).^2);
                    points.distances = distances(:);
                    points = sortrows(points,3);
                    nearestPoints = points(1:kk,:);

                    %sort by angles
                    angles = NaN(size(nearestPoints,1),1);

                    if step == 2

                        n1 = (table2array(nearestPoints(:,1:2)) - table2array(currentPoint(:,1:2))) ./ norm(table2array(nearestPoints(:,1:2)) - table2array(currentPoint(:,1:2)));  % Normalized vectors
                        n2 = (table2array(previousPoint(:,1:2)) - table2array(currentPoint(:,1:2))) / norm(table2array(previousPoint(:,1:2)) - table2array(currentPoint(:,1:2)));
                        for i = 1 : size(angles,1)
                            angles(i) = atan2(norm(det([n2; n1(i,:)])), dot(n1(i,:), n2))*180/pi;
                        end

                    else

                        n1 = [(nearestPoints.x - currentPoint.x),(nearestPoints.y - currentPoint.y)];
                        n2 = [(previousPoint.x - currentPoint.x),(previousPoint.y - currentPoint.y)];
                        for i = 1 : size(angles,1)
                            angles(i) = - atan2(det([n2; n1(i,:)]), dot(n1(i,:), n2))*180/pi;
                        end
                        angles(angles<0) = +180+(angles(angles<0)+180);

                    end

                    nearestPoints.angles = angles;
                    candidatePoints = sortrows(nearestPoints,4);%candidates

                    %verify whether all the candidates are already part of the hull
                    itsCandidate = false(size(candidatePoints,1),1);
                    for c = 1 : size(candidatePoints,1)
                        [xi,~] = polyxpoly([candidatePoints.x(c),currentPoint.x] ,[candidatePoints.y(c),currentPoint.y],hull.x,hull.y);
                        itsCandidate(c) = numel(xi)>1;
                        if ismember(candidatePoints(c,1:2), firstPoint(1,1:2))
                            itsCandidate(c) = false;
                        end
                    end

                    its = length(find(itsCandidate==true)) == size(candidatePoints,1);
                    %if at least one candidate is not--> exit from the for loop
                    if ~its
                        break
                    end

                end %for to find potential candidate points

                if ~its %whether not all the candidate points are already part of the hull
                    candidatePoints(itsCandidate==true,:) = [];

                    %find the candidate associated with the maximum angle
                    [~,indCandidate] = max(candidatePoints.angles);

                    if ~isempty(indCandidate)

                        %remove the new point from the point list
                        previousPoint = currentPoint;
                        currentPoint = candidatePoints(indCandidate,:);
                        indEquals = find(sqrt((points.x - currentPoint.x).^2 + (points.y - currentPoint.y).^2) == 0);
                        points(indEquals,:) = [];
                        %updating hull
                        hull(indHull,:) = currentPoint;
                        indHull = indHull + 1;

                        %step update
                        step = step + 1;

                    end

                else
                    %it ends
                    points = [];
                    break
                end

            else %kk > points

                %sort by angles
                angles = NaN(size(nearestPoints,1),1);

                n1 = [(nearestPoints.x - currentPoint.x),(nearestPoints.y - currentPoint.y)];
                n2 = [(previousPoint.x - currentPoint.x),(previousPoint.y - currentPoint.y)];
                for i = 1 : size(angles,1)
                    angles(i) = - atan2(det([n2; n1(i,:)]), dot(n1(i,:), n2))*180/pi;
                end
                angles(angles<0) = +180+(angles(angles<0)+180);

                nearestPoints.angles = angles;
                candidatePoints = sortrows(nearestPoints,4);%candidates

                %verify whether all the candidates are already part of the hull
                itsCandidate = false(size(candidatePoints,1),1);
                for c = 1 : size(candidatePoints,1)
                    [xi,~] = polyxpoly([candidatePoints.x(c),currentPoint.x] ,[candidatePoints.y(c),currentPoint.y],hull.x,hull.y);
                    itsCandidate(c) = numel(xi)>1;
                    if ismember(candidatePoints(c,1:2), firstPoint(1,1:2))
                        itsCandidate(c) = false;
                    end
                end

                candidatePoints(itsCandidate==true,:) = [];

                %find the candidate associated with the maximum angle
                [~,indCandidate] = max(candidatePoints.angles);

                if ~isempty(indCandidate)

                    %remove the new point from the point list
                    previousPoint = currentPoint;
                    currentPoint = candidatePoints(indCandidate,:);
                    indEquals = find(sqrt((points.x - currentPoint.x).^2 + (points.y - currentPoint.y).^2) == 0);
                    points(indEquals,:) = [];

                    %hull update
                    hull(indHull,:) = currentPoint;
                    indHull = indHull + 1;

                    %step update
                    step = step + 1;
                end

            end

        end %while points

        %check if all points inside
        allInside = false;
        if numel(find(~isnan(hull.x)))>=4
        in = inpolygon(points0.x,points0.y,hull.x,hull.y);
        [selfInt,~,~] = selfintersect(hull.x(2:end), hull.y(2:end));
        allInside = ~length(find(in==0))>0 && isempty(selfInt);
        end

    end %while all inside

    if kk0 == NPoints
        hullID = boundary(table2array(points0(:,1:2)),0.85);
        hull = points0(hullID,:);
    end

end % if points>threshold

%remove exceeding hull rows
computedPolygon = hull(find(~isnan(hull.x)),:);

end