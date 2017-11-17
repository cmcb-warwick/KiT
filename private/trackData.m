function [coords,time,coordIdx,coordsStd] = trackData(idx,dataStruct,trackStats,rotate)
%% read track coordinates

if nargin<4
  rotate=1;
end

%get indices of feature making track
featIndx = dataStruct.tracks(idx).tracksFeatIndxCG;

%get start time and end time of track
startTime = dataStruct.tracks(idx).seqOfEvents(1,1);
endTime = dataStruct.tracks(idx).seqOfEvents(2,1);
lifeTime = endTime - startTime + 1;

% read track coordinates and their stds for idx
coords = NaN(lifeTime,3);
coordsStd = NaN(lifeTime,3);
if rotate && ~isempty(dataStruct.planeFit)
    for iFrame = 1 : lifeTime
        jFrame = startTime + iFrame - 1;
        featIndxTmp = featIndx(iFrame);
        if featIndxTmp ~= 0
            coords(iFrame,1:3) = dataStruct.planeFit(jFrame).rotatedCoord(featIndxTmp,1:3);
            coordsStd(iFrame,1:3) = dataStruct.planeFit(jFrame).rotatedCoord(featIndxTmp,4:6);
        end
    end
else
    for iFrame = 1 : lifeTime
        jFrame = startTime + iFrame - 1;
        featIndxTmp = featIndx(iFrame);
        if featIndxTmp ~= 0
            coords(iFrame,1:3) = dataStruct.initCoord(jFrame).allCoord(featIndxTmp,1:3);
            coordsStd(iFrame,1:3) = dataStruct.initCoord(jFrame).allCoord(featIndxTmp,4:6);
        end
    end
end

% coordsInfo = reshape(dataStruct.tracks(idx).tracksCoordAmpCG,8,[])';
% coords = coordsInfo(:,1:3);
% coordsStd = coordsInfo(:,5:7);

% read timepoints of the track
time = (trackStats(1,1,idx):trackStats(2,1,idx))';

% remove gaps
time(all(isnan(coords),2)) = [];

% remember indices into colCoords that correspond to the timepoints
coordIdx = time - trackStats(1,1,idx) + 1;
