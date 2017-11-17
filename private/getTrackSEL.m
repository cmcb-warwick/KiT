function trackSEL = getTrackSEL(trackedFeatureInfo)
%GETTRACKSEL outputs track start times, end times and lifetimes
%
%SYNOPSIS trackSEL = getTrackSEL(trackedFeatureInfo);
%
%INPUT  trackedFeatureInfo: -- EITHER -- 
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%                           -- OR -- 
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%
%OUTPUT trackSEL          : An array with 3 columns and number of rows equal
%                           to number of (compound) tracks. 1st column 
%                           indicates track start times, 2nd column
%                           indicates track end times and 3rd column
%                           indicates track lifetimes.
%
%REMARKS The details of compound tracks are currently ignored. The start
%times, end time and life times are for the overall compound tracks, not
%their branches.
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackSEL = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getTrackSEL: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track information extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of tracks
if isstruct(trackedFeatureInfo)
    numTracks = length(trackedFeatureInfo);
else
    numTracks = size(trackedFeatureInfo,1);
end

%alocate memory for output
trackSEL = zeros(numTracks,3);

%if input is in structure format
if isstruct(trackedFeatureInfo)

    %find track start times
    for i=1:numTracks
        trackSEL(i,1) = trackedFeatureInfo(i).seqOfEvents(1,1);
    end

    %find track end times
    for i=1:numTracks
        trackSEL(i,2) = trackedFeatureInfo(i).seqOfEvents(end,1);
    end

else %if input is a matrix

    %make new matrix which contains only one column per time point
    trackedFeatureInfo = trackedFeatureInfo(:,1:8:end);
    
    %find non-empty tracks
    indxGood = find(~isnan(max(trackedFeatureInfo,[],2)));

    %find track start times
    for i = indxGood'
        trackSEL(i,1) = find(~isnan(trackedFeatureInfo(i,:)),1,'first');
    end

    %find track end times
    for i = indxGood'
        trackSEL(i,2) = find(~isnan(trackedFeatureInfo(i,:)),1,'last');
    end

end %(if isstruct(trackedFeatureInfo))
    
%calculate track lifetimes
trackSEL(:,3) = trackSEL(:,2) - trackSEL(:,1) + 1;


%%%%% ~~ the end ~~ %%%%%

