function kitShowDragonTails(job,varargin)
% KITSHOWDRAGONTAILS
%

% default options
opts.channel = 1;
opts.contrast = {[0.1 1],[0.1 1],[0.1 1]};
opts.tailLength = 10; %timepoints
opts.timePoint = 1; %centre of dragon tail
opts.tracks = [];
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

% get movie length
movieLength = job.metadata.nFrames;
if opts.timePoint > movieLength
  error('Please provide a central timepoint within the length of the movie (max %i time points).',movieLength);
end

% create time series for plotting
t = max(opts.timePoint-opts.tailLength,1): ...
    min(opts.timePoint+opts.tailLength,movieLength);

% check whether this track is within this time frame for greater than 2 time points
for iTrack = opts.tracks
  trackT = tracks(iTrack).seqOfEvents(1,1):tracks(iTrack).seqOfEvents(2,1);
  if length(intersect(trackT,t)) < 3
    opts.tracks = setdiff(opts.tracks,iTrack);
  end
end

% show the image
kitShowImage(job,'imageChans',opts.channel,'timePoint',opts.timePoint,'contrast',opts.contrast)
hold on

for iTrack = opts.tracks
  
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
  plot(trackCoords(1:opts.tailLength+1,1),trackCoords(1:opts.tailLength+1,2),'g')
  plot(trackCoords(opts.tailLength+1:end,1),trackCoords(opts.tailLength+1:end,2),'r')


end