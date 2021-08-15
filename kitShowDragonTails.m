function kitShowDragonTails(job,varargin)
% KITSHOWDRAGONTAILS To visualise tracks on a still image
%
% KITSHOWDRAGONTAILS(job, ...) Plots tracked kinetochores on top of a movie still
% 
% Options, defaults in {}:-
% 
%    channel: {1}, 2 or 3. Which channel to show.
%
%    contrast: {[0.1 1]} or similar two-element vector. Range over which to
%           contrast images. Tips:
%               - Increase brightness by changing to [0.1 0.9]
%               - Decrease background noise by changing to [0.5 1]
%
%    crop: {1}, 0 or -1. Whether or not to show cropped images, as defined
%       by ROI.crop. -1 shows a full image with the cropped region
%       annotated. Passed to kitShowImage. 
%    
%    tailLength: {10} or number. Length of tracks to plot prior to and after selected frame.
%
%    timePoint: {1} or number. Timepoint at which to plot spots.
%
%    tracks: {[]} or array of numbers. The indices of tracks or pairs to plot from a movie.
%           These numbers will refer to the track ID for each track.
%    
%    usePairs: {0}, 1. Whether to only plot tracks for paired sister kinetochores,
%           or for all tracked kinetochores. Pairs will be indexed by Sister ID 
%            rather than track ID.
%
% Example usage: kitShowDragonTails(job,'crop',-1,'tracks',[58,93,133,143,168,197],'tailLength',99,'timePoint',1)
%
% Copyright (c) 2020 Jonathan U Harrison
%%%%%%%%%%%%%%%%%%%%%%%%%%

% default options
opts.channel = 1;
opts.contrast = {[0.1 1],[0.1 1],[0.1 1]};
opts.crop = 1;
opts.plotColors = {'c','m'};
opts.tailLength = 10; %timepoints
opts.timePoint = 1; %centre of dragon tail
opts.tracks = [];
opts.usePairs = 0;
% process user-defined options
opts = processOptions(opts,varargin{:});

% get important data
dataStruct = job.dataStruct{opts.channel};
initCoord = dataStruct.initCoord;
tracks = dataStruct.tracks;

% check which tracks to plot - ensure no absent tracks have been requested
if isempty(opts.tracks)
    opts.tracks = 1:length(tracks);
else
    opts.tracks = intersect(opts.tracks,1:length(tracks));
end

trackPairs = dataStruct.sisterList(1).trackPairs(:,1:2);
%subset = trackPairs(opts.tracks,1:2); subset = subset(:)';

% get movie length
movieLength = job.metadata.nFrames;
if opts.timePoint > movieLength
    error('Please provide a central timepoint within the length of the movie (max %i time points).',movieLength);
end

% create time series for plotting
t = max(opts.timePoint-opts.tailLength,1): ...
    min(opts.timePoint+opts.tailLength,movieLength);

% show the image
kitShowImage(job,'imageChans',opts.channel,'timePoint',opts.timePoint,'contrast',opts.contrast,'crop',opts.crop)
hold on
if opts.usePairs
    for iPair = opts.tracks
        trackPairs = dataStruct.sisterList(1).trackPairs(:,1:2);
        plot_track_on_movie(t,tracks,initCoord,opts,trackPairs(iPair,1));
        plot_track_on_movie(t,tracks,initCoord,opts,trackPairs(iPair,2));
    end
else
    % check whether this track is within this time frame for greater than 2 time points
    for iTrack = opts.tracks
        trackT = tracks(iTrack).seqOfEvents(1,1):tracks(iTrack).seqOfEvents(2,1);
        if length(intersect(trackT,t)) < 3
            opts.tracks = setdiff(opts.tracks,iTrack);
        end
    end
    for iTrack = opts.tracks
        plot_track_on_movie(t,tracks,initCoord,opts,iTrack);
    end
end
end

    function plot_track_on_movie(t,tracks,initCoord,opts,iTrack)
    % check whether this track is within this time frame for greater than 2 time points
    trackT = tracks(iTrack).seqOfEvents(1,1):tracks(iTrack).seqOfEvents(2,1);
    if length(intersect(trackT,t)) < 3
        opts.tracks = setdiff(opts.tracks,iTrack);
    end
        % predesignate trackCoords
        trackCoords = nan(length(t),3);
        featIdxs = nan(length(t),1);
        
        % get timing information
        start = tracks(iTrack).seqOfEvents(1,1);
        fin   = tracks(iTrack).seqOfEvents(2,1);
        trackT = max(t(1),start):min(t(end),fin);
        vectT = trackT-t(1)+1;
        trackT = trackT-start+1;
        % and featIdxs
        featIdxs(vectT) = tracks(iTrack).tracksFeatIndxCG(trackT);
        
        % get pixel coordinate information
        for iTime = 1:length(t)
            if ~isnan(featIdxs(iTime)) && featIdxs(iTime)>0
                trackCoords(iTime,:) = initCoord(t(iTime)).allCoordPix(featIdxs(iTime),1:3);
            end
        end
        
        % plot tracks
        plot(trackCoords(1:opts.tailLength+1,1),trackCoords(1:opts.tailLength+1,2),opts.plotColors{1},'LineWidth',1.5)
        plot(trackCoords(opts.tailLength+1:end,1),trackCoords(opts.tailLength+1:end,2),opts.plotColors{2},'LineWidth',1.5)
    end
