function dataStruct = kitGroupSisters(dataStruct,verbose,allowMulti)
% KITGROUPSISTERS groups sister kinetochores
%
% SYNOPSIS: dataStruct = kitGroupSisters(dataStruct)
%
% INPUT dataStruct: data structure as created by kitMakeDataStruct, with
%                   the fields: "dataProperties", "initCoord", "planeFit" &
%                   "tracks". Field "planeFit" can be empty.
%       verbose (opt) : 0 - no plotting (default)
%                       1 - plot 4 frames with sister assignment
%                       2 - 1 & plot all tracks
%       allowMulti: boolean, allow multiple connections between grouped
%       sisters. Optimization performed via mosek library. 
%
% OUTPUT dataStruct: Same as input, but with field sisterList
%              sisterList is a structure with length equal to the number of
%              sister kinetochore pairs.
%              sisterList(1).trackPairs is an nPairs-by-6 array with
%                   [track1,track2,cost,avg. dist,variance,alignment], that
%                   is sorted according to increasing cost.
%                   track1,2: track indices as in dataStruct.tracks.
%                   cost: cost of grouping
%                   avg. dist: average distance between the two tracks
%                   variance: variance of the distance between the tracks
%                   alignment: f(tan(alpha)), where alpha is the average
%                       angle between the distanceVector and the first
%                       eigenVector of planeFit.eigenVectors
%              sisterList(iPair).coords1 is a nTimepoints-by-6 array with
%                   the coordinates of the first of the two tracks and its
%                   std.
%              sisterList(iPair).coords2 is a nTimepoints-by-6 array with
%                   the coordinates of the second of the two tracks and its
%                   std.
%              sisterList(iPair).sisterVectors is a nTimepoints-by-6 
%                   array with the vector connecting the two sisters and
%                   its std.
%              sisterList(iPair).distances is a nTimepoints-by-2 array with
%                   the distance between sisters and its std.
%
% REMARKS Sister identification is based on globally minimizing (1) the
%           average distance between sisters, (2) the variance of the distance
%           between sisters, and (3) the alignment of sisters with the normal to
%           the metaphase plate (if relevant). 
%         Anaphase frames are not used in sister identification.
%         At the end of the code is a plotting function for the distance
%           between the tracks for debugging
%         The code cannot handle merged/splitted tracks!
%
% Copyright (c) 2007 Jonas Dorn, Khuloud Jaqaman
% Edited Jonathan U Harrison 2019-03-20

%% TEST INPUT & READ PARAMETERS

if nargin < 2 || isempty(verbose)
    verbose = 0;
end

if nargin <3
    allowMulti=0;
end

% read movieLength
nTimepoints = dataStruct.dataProperties.movieSize(4);

% read parameters
minOverlap = min(dataStruct.dataProperties.groupSisters.minOverlap,round(nTimepoints/10));
% if minOverlap < 10
%     minOverlap = 10;
% end
maxDist = dataStruct.dataProperties.groupSisters.maxDist;
maxAngle = dataStruct.dataProperties.groupSisters.maxAngle * pi / 180;
robust = dataStruct.dataProperties.groupSisters.robust;
useAlignment = dataStruct.dataProperties.groupSisters.useAlignment;
useAnaphase = dataStruct.dataProperties.groupSisters.useAnaphase;

% read track statistics. This will work only if no merge/split, i.e. if
% there are only two events per track: a start and a finish
try
    trackStats = catStruct(3,'dataStruct.tracks.seqOfEvents');
catch
    error('kitGroupSisters cannot handle merging/splitting')
end

% get track lengths
trackLength = squeeze(trackStats(2,1,:)-trackStats(1,1,:)+1);

% select tracks whose length is larger than the minimum overlap
goodTracks = find(trackLength>=minOverlap);
nGoodTracks = length(goodTracks);

%% READ TRACK INFORMATION

% preassign matrices
[variances,distances,alignment] = deal(NaN(nGoodTracks));

%find frames that have a plane (to calculate alignment cost)
%if none of the frames have a plane, we cannot use the alignment criterion
framesWiPlane = [];
for t = 1 : nTimepoints
    if ~isempty(dataStruct.planeFit) && ~isempty(dataStruct.planeFit(t).planeVectors)
        framesWiPlane = [framesWiPlane; t];
    end
end
if isempty(framesWiPlane)
    useAlignment = 0;
end

%find anaphase frames
if ~isempty(dataStruct.planeFit)
    framePhase = vertcat(dataStruct.planeFit.phase);
else
    framePhase = repmat('e',nTimepoints,1);
end
anaphaseFrames = find(framePhase == 'a');
if isempty(anaphaseFrames)
    lastFrameNotAna = nTimepoints;
else
    lastFrameNotAna = anaphaseFrames(1) - 1;
end
%if the whole movie is in anaphase, there's no point in looking for
%sisters. Exit with an empty sisterList.
if length(anaphaseFrames) == nTimepoints
    sisterList = struct('trackPairs',[],'coords1',[],...
        'coords2',[],'sisterVectors',[],'distances',[]);
    dataStruct.sisterList = sisterList;
    return
end

% read normals to plane
normals = NaN(nTimepoints,3);
if useAlignment == 1 && ~isempty(dataStruct.planeFit)
    normals(framesWiPlane,:) = catStruct(2,'dataStruct.planeFit.planeVectors(:,1)')';
end

% loop through the good tracks, calculate for every pair mean distance
% and variance
if verbose == 2
    figure
    hold on
end
for jTrack = 1:nGoodTracks % loop cols

    % read index of track
    jIdx = goodTracks(jTrack);

    % read track coordinates etc. (one could get the coordinate stds as
    % well, but for now no need here - KJ)
    % coordinates are in metaphase plane rotated frame of reference
    [colCoords,colTime,colIdx] = trackData(jIdx,dataStruct,trackStats);

    % plot individual tracks
    if verbose == 2
        plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),...
            'Color',extendedColors(jTrack))
    end

    for iTrack = jTrack+1:nGoodTracks % loop rows

        % read index of track
        iIdx = goodTracks(iTrack);

        % read track coordinates (one could get the coordinate stds as
        % well, but for now no need - KJ)
        % coordinates are in metaphase plane rotated frame of reference
        [rowCoords,rowTime,rowIdx] = trackData(iIdx,dataStruct,trackStats);

        % find common time
        [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

        %separate common time which is in anaphase from common time which
        %is not in anaphase
        ctNotAna = find(framePhase(commonTime)~='a');
        ctAna = find(framePhase(commonTime)=='a');
        commonTimeNA = commonTime(ctNotAna); %not in anaphase
        ctColIdxNA = ctColIdx(ctNotAna);
        ctRowIdxNA = ctRowIdx(ctNotAna);
        commonTimeA = commonTime(ctAna); %in anaphase
        ctColIdxA = ctColIdx(ctAna);
        ctRowIdxA = ctRowIdx(ctAna);
        clear commonTime ctColIdx ctRowIdx
        
        %if the common time between the two tracks is at least minOverlap,
        %calculate parameters (otherwise, they stay as NaN, which assigns 
        %them -1 in linkTracks)
        if length(commonTimeNA) >= minOverlap
            
            %get the positions of the two tracks at the end of anaphase
            %and add up their signs
            if ~isempty(ctAna)
                track1Pos = colCoords(colIdx(ctColIdxA(end)),1);
                track2Pos = rowCoords(rowIdx(ctRowIdxA(end)),1);
                sumSigns = sign(track1Pos) + sign(track2Pos);
            else
                sumSigns = 0;
            end
            
            %if there is anaphase, try to link pairs only if they end up on
            %opposite sides of the metaphase plate (if requested by user)
            %(otherwise, they stay as NaN, which assigns them -1 in linkTracks)
            if (useAnaphase && (sumSigns==0)) || ~useAnaphase

                % calculate distance (microns)
                distanceVector = colCoords(colIdx(ctColIdxNA),:) -...
                    rowCoords(rowIdx(ctRowIdxNA),:);

                % calculate alignment:

%                 % no need for this now since the coordinates are in the
%                 % metaphase-plate frame of reference
%                 % retain commonTime-normals
%                 commonNormals = normals(commonTimeNA,:);
% 
%                 % cos(alpha) = dot(a,b)/(norm(a)*norm(b))
%                 [distance,distanceVectorN] = normList(distanceVector);
%                 alpha = acos(abs(dot(distanceVectorN',commonNormals')));

                %get the angle between distance vector and normal
                [distance,distanceVectorN] = normList(distanceVector);
                alpha = acos(abs(distanceVectorN(:,1)));

                % average alpha, rather than tan to be nice to pairs that will
                % align eventually. Potentially this can be put under the control
                % of the "robust" switch, too
                %average alpha only over frames where there is a plane (the
                %rest are NaN). If none of the frames have a plane, the average
                %will be NaN.
                %also get the standard deviation of alpha
                meanAlpha = nanmean(alpha);
                stdAlpha = nanstd(alpha);

                % get distance mean and standard deviation
                if robust
                    [rMean,rStd]=robustMean(distance);
                else
                    rMean = mean(distance);
                    rStd = std(distance);
                end
                %assign distance mean for pair
                distances(iTrack,jTrack) = rMean;
                distances(jTrack,iTrack) = rMean;
                %assign distance variance for pair
                variances(jTrack,iTrack) = rStd^2;
                variances(iTrack,jTrack) = rStd^2;

                %assign alignment cost for pair if the average angle is less
                %than maxAngle degrees. Otherwise, keep as NaN to prohibit the link
                if meanAlpha < maxAngle
		            if verbose
                        fprintf('\nrobust mean distance is %f pm %f compared to %f\n',rMean,rStd,maxDist);
                        fprintf('\nrobust mean allignment is %f pm %f compared to %f\n',meanAlpha,stdAlpha,maxAngle);
                    end
                    [alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
                        = deal(2*sqrt(3)*tan(meanAlpha)+1);
                end
                
            end %(if (useAnaphase && (sign(track1Pos)+sign(track2Pos))==0) || ~useAnaphase)

        end %(if length(commonTime) >= minOverlap)
        
    end %(for iTrack = jTrack+1:nGoodTracks)
end %(for jTrack = 1:nGoodTracks)

%% CREATE COST MATRIX & GROUP
[r2c,c2r,costMat,linkedIdx] = ...
    linkTracks(distances,variances,alignment,...
    nGoodTracks,maxDist,useAlignment,allowMulti,verbose);

if all(isnan(r2c))
    sisterList = struct('trackPairs',[],'coords1',[],...
        'coords2',[],'sisterVectors',[],'distances',[]);
    dataStruct.sisterList = sisterList;
    return;
end

%find pairs that get a unique assignment
if allowMulti
    c2rmosek = int32(c2r);
    c2r=r2c;
end
goodPairIdxL = r2c==c2r;

if verbose
    % plot for 4 frames (that's about how many can be properly displayed)
    deltaT = max(1,floor(nTimepoints/4));
    tOffset = max(1,ceil(deltaT/2));
    t=tOffset:deltaT:nTimepoints;
    % highlight polygons
    r2cTmp = r2c;
    r2cTmp(~goodPairIdxL) = -r2cTmp(~goodPairIdxL);
    plotGroupResults(t,r2cTmp,nGoodTracks,...
        goodTracks,dataStruct,distances,maxDist,...
        'Initial grouping for %s. G/B-Cutoff=distance');
end

%% RESOLVE POLYGONS

% identify polygons. Polygons have to be closed, thus, it should not matter
% where we start. Also, since the distance between polygons and the rest is
% hopefully fairly large, we don't care about neighborhood.
polygonIdx = find(~goodPairIdxL);

% remove the not-linked tracks
polygonIdx(isnan(r2c(polygonIdx))) = [];
polyList = [];
while ~isempty(polygonIdx)
    polyList(1) = polygonIdx(1);
    polygonIdx(1) = [];
    done = false;
    while ~done
        % look up the row the last corner links to
        nextCorner = r2c(polyList(end));
        % check whether the new corner has already been used
        if any(nextCorner == polyList)
            % if yes, exit. The polygon is complete
            done = true;
        else
            polyList(end+1) = nextCorner; %#ok<AGROW>
            % remove corner from polygonIdx
            polygonIdx(polygonIdx==nextCorner) = [];
        end
    end % identify polygon

    % within the polygon: find closest distance to identify first pair.
    % Remove it, and check for more pairs. This will potentially result in
    % more pairs than removal of large distances starting from a tetragon
    done = false;
    while ~done
        % read current cost matrix
        currentCost = costMat(polyList,polyList);
        currentCost(currentCost==-1) = inf;
        % find pair with lowest cost
        [v,minIdx] = min(currentCost(:));
        [idx1,idx2] = ind2sub(size(currentCost),minIdx);
        % write pair into r2c
        r2c(polyList(idx1)) = polyList(idx2);
        r2c(polyList(idx2)) = polyList(idx1);
        c2r(polyList(idx1)) = polyList(idx2);
        c2r(polyList(idx2)) = polyList(idx1);

        % check whether there are still tracks to link
        polyList([idx1,idx2]) = [];

        if length(polyList) > 1 && any(~isinf(currentCost(:)))
            % continue
        else
            % clear polyList, write NaN into r2c, c2r
            r2c(polyList) = NaN;
            c2r(polyList) = NaN;
            polyList = [];
            done = true;
        end
    end % resolve individual polygons

end % resolve all polygons

if verbose
    % plot final version
    r2cTmp = r2c;
    r2cTmp(~goodPairIdxL) = -r2cTmp(~goodPairIdxL);
    plotGroupResults(t,r2cTmp,nGoodTracks,...
        goodTracks,dataStruct,costMat,maxDist,...
        'Final grouping for G/B-Cutoff=cost');
end



%% assemble sister information

% sisterList:
%   .trackPairs
%   .coords1
%   .coords2
%   .sisterVectors
%   .distances

if allowMulti
    nGoodPairs = sum(c2rmosek>0 & ~isnan(r2c));
sisterList(1:nGoodPairs/2,1) = ...
    struct('trackPairs',[],'coords1',NaN(nTimepoints,6),...
    'coords2',NaN(nTimepoints,6),'sisterVectors',NaN(nTimepoints,6),...
    'distances',NaN(nTimepoints,2));

% write trackPairs. Store: pair1,pair2,cost,dist,var,alignment
goodPairIdxL = (c2rmosek>0) & ~isnan(r2c);
sisterList(1).trackPairs = ...
    [goodTracks(c2rmosek(goodPairIdxL)),goodTracks(r2c(goodPairIdxL)),...
    costMat(linkedIdx),distances(linkedIdx),variances(linkedIdx),...
    alignment(linkedIdx)];
sisterList(1).trackPairs(isnan(sisterList(1).trackPairs(:,6)),6) = 0;
else
    goodPairIdxL = r2c==c2r;
    linkedIdx=sub2ind([nGoodTracks nGoodTracks],find(goodPairIdxL),r2c(goodPairIdxL));
    nGoodPairs = sum(goodPairIdxL);
sisterList(1:nGoodPairs/2,1) = ...
    struct('trackPairs',[],'coords1',NaN(nTimepoints,6),...
    'coords2',NaN(nTimepoints,6),'sisterVectors',NaN(nTimepoints,6),...
    'distances',NaN(nTimepoints,2));

% write trackPairs. Store: pair1,pair2,cost,dist,var,alignment
    sisterList(1).trackPairs = ...
    [goodTracks(goodPairIdxL),goodTracks(r2c(goodPairIdxL)),...
    costMat(linkedIdx),distances(linkedIdx),variances(linkedIdx),...
    alignment(linkedIdx)];
sisterList(1).trackPairs(isnan(sisterList(1).trackPairs(:,6)),6) = 0;
end

% remove redundancy
sisterList(1).trackPairs(:,1:2) = sort(sisterList(1).trackPairs(:,1:2),2);
sisterList(1).trackPairs = unique(sisterList(1).trackPairs,'rows');

% sort according to cost
sisterList(1).trackPairs = sortrows(sisterList(1).trackPairs,3);

%need to combine groups together into single tracks
nSisterPairs = size(sisterList(1).trackPairs,1); %counting duplicates to merge
[~,ISister1,~] = unique(sisterList(1).trackPairs(:,1), 'rows');
ixDupSister1 = setdiff(1:nSisterPairs, ISister1);
dupSister1Values = sisterList(1).trackPairs(ixDupSister1,1);

[~,ISister2,~] = unique(sisterList(1).trackPairs(:,2), 'rows');
ixDupSister2 = setdiff(1:nSisterPairs, ISister2);
dupSister2Values = sisterList(1).trackPairs(ixDupSister2,2);


% loop over trackPairs to get their coordinates and distances
for iPair = 1:nSisterPairs 
    %%%%%%%%%%%
    [lia1,lib1] = ismember(iPair,ixDupSister1); %check for duplicates in first column
    [lia2,lib2] = ismember(iPair,ixDupSister2); %check for duplicates in second column
    if lia1
    %locate duplicates
    dupPairIdx = find(sisterList(1).trackPairs(:,1)==dupSister1Values(lib1));
    %get information for new second sister
    [colCoords,colTime,colIdx,colCoordsStd] = ...
        trackData(sisterList(1).trackPairs(iPair,2),dataStruct,trackStats);
    sisterList(dupPairIdx(1)).coords2(colTime,:) = [colCoords(colIdx,:) ...
        colCoordsStd(colIdx,:)];
    elseif lia2 
            %locate duplicates
    dupPairIdx = find(sisterList(1).trackPairs(:,2)==dupSister2Values(lib2));
    %get information for new first sister
    [rowCoords,rowTime,rowIdx,rowCoordsStd] = ...
        trackData(sisterList(1).trackPairs(iPair,1),dataStruct,trackStats);
    sisterList(dupPairIdx(1)).coords1(rowTime,:) = [rowCoords(rowIdx,:) ...
        rowCoordsStd(rowIdx,:)];
    else %is not repeat occurrence of a duplicate  
    %%%%%%%%%%
    
    %get information for first sister
    [rowCoords,rowTime,rowIdx,rowCoordsStd] = ...
        trackData(sisterList(1).trackPairs(iPair,1),dataStruct,trackStats);
    
    %get information for second sister
    [colCoords,colTime,colIdx,colCoordsStd] = ...
        trackData(sisterList(1).trackPairs(iPair,2),dataStruct,trackStats);

    %find common time between them
    [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

%EDITED SO THAT WE DO NOT THROW AWAY INFORMATION ABOUT THE LONGER TRACK -JUH
    %store the coordinates of the first sister
    sisterList(iPair).coords1(rowTime,:) = [rowCoords(rowIdx,:) ...
        rowCoordsStd(rowIdx,:)];

    %store the coordinates of the second sister
    sisterList(iPair).coords2(colTime,:) = [colCoords(colIdx,:) ...
        colCoordsStd(colIdx,:)];

    %calculate the vector connecting the two sisters and its std (microns)
    sisterVectors = [colCoords(colIdx(ctColIdx),:) - rowCoords(rowIdx(ctRowIdx),:) ...
        sqrt(colCoordsStd(colIdx(ctColIdx),:).^2 + rowCoordsStd(rowIdx(ctRowIdx),:).^2)];
    sisterList(iPair).sisterVectors(commonTime,:) = sisterVectors;
    
    %calculate the distance between the two sisters and its std (microns)
    sisterDist = sqrt(sum(sisterVectors(:,1:3).^2,2));
    sisterDistStd = sqrt(sum((sisterVectors(:,1:3)./repmat(sisterDist,1,3)).^2 .* ...
        sisterVectors(:,4:6).^2,2));
    sisterList(iPair).distances(commonTime,:) = [sisterDist sisterDistStd];
    end
end % loop goodPairs
sisterList([ixDupSister1,ixDupSister2])=[];
nGoodPairs = nGoodPairs - 2*length([ixDupSister1,ixDupSister2]);
%% remove extra large distances from sister pairing

%for math behind algorithm, see Danuser, 1992 or Rousseeuw & Leroy, 1987

%put all sister distances in one vector
sisterDist = [];
for iSister = 1 : nGoodPairs/2
    sisterDist = [sisterDist; sisterList(iSister).distances(1:lastFrameNotAna,1)];
end
sisterDist = sisterDist(~isnan(sisterDist));

%calculate median of distances
medDistance = median(sisterDist);

%get residuals, i.e. distance of observations from median
residuals = sisterDist - medDistance;

%square the residuals
res2 = residuals .^ 2;

%calculate the median of the squared residuals
medRes2 = median(res2);

%assign a res2 of zero to all distances smaller than the median because we
%do not want to remove any of those
res2(residuals < 0) = 0;

%define parameters to remove outliers
k = 2.5; %assuming Gaussian distribution, keeps 95% of the distances
magicNumber2 = 1.4826^2; %see same publications

%calculate testvalues to determine which observations are inliers and which
%are outliers
%distances smaller than median will automatically get a testValue of zero
%so that they would be considered inliers
testValue = res2 / (magicNumber2 * medRes2);

%calculate the largest inlier distance
maxSisterDist = max(sisterDist(testValue <= k^2));

%remove all distances larger than this maximum distance
for iPair = 1 : nGoodPairs/2
    
    %find indices of distances larger than maximum allowed
    outlierIndx = find(sisterList(iPair).distances(1:lastFrameNotAna,1) ...
        > maxSisterDist);
    
    %remove those timepoints from the sister information
    sisterList(iPair).coords1(outlierIndx,:) = NaN;
    sisterList(iPair).coords2(outlierIndx,:) = NaN;
    sisterList(iPair).sisterVectors(outlierIndx,:) = NaN;
    sisterList(iPair).distances(outlierIndx,:) = NaN;
    
end %(for iPair = 1 : nGoodPairs/2)

%% assign output to dataStruct

dataStruct.sisterList = sisterList;


%% SUBFUNCTIONS & DEBUG HELPER FUNCTIONS

%% link tracks
function [r2c,c2r,costMat,linkedIdx] = linkTracks(distances,variances,...
    alignment,nGoodTracks,maxDist,useAlignment,allowMulti, verbose)

if verbose
% cutoff distances
figure; 
[f,x] = ecdf(distances(:));
stairs(x,f,'linewidth',3);
title('Distances');
set(gca,'fontsize',20);
figure;
[f,x] = ecdf(variances(:));
stairs(x,f,'linewidth',3);
title('Variances');
set(gca,'fontsize',20);
figure;
[f,x] = ecdf(alignment(:));
stairs(x,f,'linewidth',3);
title('Alignment');
set(gca,'fontsize',20);
end
distCutoffIdx = distances>maxDist;
distances(distCutoffIdx) = NaN;
variances(distCutoffIdx) = NaN;
alignment(distCutoffIdx) = NaN;

% make cost matrix
switch useAlignment
    case 0
        costMat = distances.*variances;
    case 1
        costMat = distances.*variances.*alignment;
end

% replace NaNs with -1
costMat(isnan(costMat)) = -1;
sum(costMat(:)~=-1)
% lap costMat
if all(costMat==-1)
    r2c = NaN(2*nGoodTracks,1);
    c2r = NaN(2*nGoodTracks,1);
else
    if allowMulti
        tic;    [r2c,c2r] = lapMosek(costMat,-1,1,1); toc
        fprintf('%d extra connections included \n',sum(1 - diff(c2r)));
        %c2r = r2c;        
        %tic;    [r2c,c2r] = lap(costMat,-1,0,1,1); toc
        %allowMulti=0;
    else
        tic;    [r2c,c2r] = lap(costMat,-1,0,1,1); toc
    end
end

if allowMulti
    % shorten r2c, c2r. No link is NaN
r2c = double(r2c);
r2c(r2c>nGoodTracks) = NaN;
c2r = double(c2r);
c2r(c2r>nGoodTracks) = NaN;
linkedIdx=sub2ind([nGoodTracks,nGoodTracks],c2r,r2c);
linkedIdx(isnan(linkedIdx)) = [];

else
    % shorten r2c, c2r. No link is NaN
r2c = double(r2c(1:nGoodTracks));
r2c(r2c>nGoodTracks) = NaN;
c2r = double(c2r(1:nGoodTracks));
c2r(c2r>nGoodTracks) = NaN;

linkedIdx=sub2ind([nGoodTracks nGoodTracks],1:nGoodTracks,r2c(1:nGoodTracks)');
linkedIdx(isnan(linkedIdx)) = [];
end

%% PLOT/DEBUG FUNCTIONS
function plotGroupResults(tt,r2c,nGoodTracks,goodTracks,dataStruct,variances,cutoff,figureName)
% plot sister-links for one frame. Number indicates goodTrackIdx
% green: variance below cutoff
% blue : variance above cutoff
% red: track with no partner

if nargin < 7 || isempty(cutoff)
    cutoff = inf;
end
figure('Name',figureName)
ntt = length(tt);
rows = ceil(sqrt(ntt));
cols = ceil(ntt/rows);
for p=1:ntt
    t = tt(p);
    [c2pg,c2pb,c2pr]=deal(nan(3*nGoodTracks,3));

    for i=nGoodTracks:-1:1
        idx1 = goodTracks(i);
        % check out of bounds
        if t<dataStruct.tracks(idx1).seqOfEvents(1,1) || t>dataStruct.tracks(idx1).seqOfEvents(2,1)
            c1 = nan(1,3);
        else
            % read via idx, not time
            tIdx = t-dataStruct.tracks(idx1).seqOfEvents(1,1)+1;
            c1=dataStruct.tracks(idx1).tracksCoordAmpCG((tIdx-1)*8+1:(tIdx-1)*8+3);
        end
        if isnan(r2c(i)) || abs(r2c(i)) > nGoodTracks
            c2=nan(1,3);
            v=NaN;
        else
            idx2 = goodTracks(abs(r2c(i)));
            if t<dataStruct.tracks(idx2).seqOfEvents(1,1) || t>dataStruct.tracks(idx2).seqOfEvents(2,1)
                c2 = nan(1,3);
            else
                % read via idx, not time
                tIdx = t-dataStruct.tracks(idx2).seqOfEvents(1,1)+1;
                c2=dataStruct.tracks(idx2).tracksCoordAmpCG((tIdx-1)*8+1:(tIdx-1)*8+3);
            end
            if r2c(i)<0
                v = nan; % polygons are also red
            else
                v = variances(i,r2c(i));
            end
        end

        if v < cutoff
            c2pg(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        elseif isnan(v)
            c2pr(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        else
            c2pb(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        end


    end
    subplot(rows,cols,p)
    plot3(c2pg(:,1),c2pg(:,2),c2pg(:,3),'*-g')
    hold on
    plot3(c2pb(:,1),c2pb(:,2),c2pb(:,3),'*-b')
    plot3(c2pr(:,1),c2pr(:,2),c2pr(:,3),'*-r')
    text(c2pg(1:3:end,1),c2pg(1:3:end,2),c2pg(1:3:end,3),num2str((1:nGoodTracks)'))
    text(c2pb(1:3:end,1),c2pb(1:3:end,2),c2pb(1:3:end,3),num2str((1:nGoodTracks)'))
    text(c2pr(1:3:end,1),c2pr(1:3:end,2),c2pr(1:3:end,3),num2str((1:nGoodTracks)'))
    grid on
    title(sprintf('Frame %i',t))
end

%% PLOT DISTANCE BETWEEN TRACKS
function distance = plotTrackDistance(iIdxT,jIdxT,goodTracks,dataStruct,trackStats) %#ok<DEFNU>
% plots distance between two tracks vs. time; two connected tracks

jIdx = goodTracks(jIdxT);
iIdx = goodTracks(iIdxT);

[rowCoords,rowTime,rowIdx] = trackData(iIdx,dataStruct,trackStats);
[colCoords,colTime,colIdx] = trackData(jIdx,dataStruct,trackStats);
% find common time
[commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

% calculate distance (microns)
distance = sqrt(sum((colCoords(colIdx(ctColIdx),:) -...
    rowCoords(rowIdx(ctRowIdx),:)).^2,2));

figure('Name',sprintf('tracks %i (r) & %i (b)',iIdxT,jIdxT))

subplot(1,2,1)
plot(commonTime,distance)
ylim([0,max(5,max(distance))])

% also plot 3d track
subplot(1,2,2)
plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),'b')
hold on
plot3(rowCoords(:,1),rowCoords(:,2),rowCoords(:,3),'r')
cc = reshape([colCoords(colIdx(ctColIdx),:),...
    rowCoords(rowIdx(ctRowIdx),:),nan(length(distance),3)]',3,[])';
plot3(cc(:,1),cc(:,2),cc(:,3),'g')
plot3(cc(1:2,1),cc(1:2,2),cc(1:2,3),'*')
