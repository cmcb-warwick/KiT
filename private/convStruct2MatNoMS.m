function [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal)
%CONVSTRUCT2MATNOMS converts tracks from structure format to matrix format, provided there are NO merges/splits.
%
%SYNPOSIS [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal)
%
%INPUT  tracksFinal: Output of trackCloseGapsKalman, when run with
%                    gapCloseParam.mergeSplit = 0.
%       trackedFeatureInfo, trackedFeatureIndx: Output of trackWithGapClosing.
%
%Khuloud Jaqaman, February 2008

%% conversion

%get number of tracks
numTracks = length(tracksFinal);

%get number of time points
tmp = vertcat(tracksFinal.seqOfEvents);
numTimePoints = max(tmp(:,1));

%reserve memory for matrix of tracks
trackedFeatureInfo = NaN(numTracks,8*numTimePoints);

%reserve memory for matrix of feature indices
trackedFeatureIndx = zeros(numTracks,numTimePoints);

%put tracks in matrix
for iTrack = 1 : numTracks
    startTime = tracksFinal(iTrack).seqOfEvents(1,1);
    endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
    trackedFeatureInfo(iTrack,8*(startTime-1)+1:8*endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG;
    trackedFeatureIndx(iTrack,startTime:endTime) = ...
        tracksFinal(iTrack).tracksFeatIndxCG;
end

%% ~~~ the end ~~~