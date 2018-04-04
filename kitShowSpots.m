function kitShowSpots(job,varargin)
% KITSHOWSPOTS Plots images of each spot in a movie.
%
%    KITSHOWSPOTS(JOB,...) Plots coordinates over images of each
%    spot localised in a given channel at a given timepoint.
%
%    Options, defaults in {}:-
%
%    channel: {1}, 2 or 3. Which channel to show.
%
%    contrast: {[0.1 1]} or similar two-element vector. Range over which to
%           contrast images. Tips:
%               - Increase brightness by changing to [0.1 0.9]
%               - Decrease background noise by changing to [0.5 1]
%
%    plotCoords: 0 or {1}. Whether or not to plot coordinates over images.
%
%    subset: {[]} or array of numbers. The spots to plot from a movie.
%           These numbers will refer to the track ID for each spot.
%
%    timePoint: {1} or number. Timepoint at which to plot spots.
%
%    zProject: {0}, 1 or -1. Whether or not to project in the z-direction.
%           -1 will project in the 5 z-slices surrounding the spot.
%
% Copyright (c) 2017 C. A. Smith

opts.channel = 1;
opts.contrast = [0.1 1];
opts.plotCoords = 1;
opts.subset = [];
opts.timePoint = 1;
opts.zProject = 0;
opts = processOptions(opts,varargin{:});

% suppress warnings
w = warning;
warning('off','all');

% get important information
dataStruct = job.dataStruct{opts.channel};
if isempty(opts.subset)
  if isfield(dataStruct,'tracks')
    nSpots = length(dataStruct.tracks);
    opts.subset = 1:nSpots;
  else
    nSpots = size(dataStruct.initCoord(opts.timePoint).allCoord,1);
    opts.subset = 1:nSpots;
  end
else
  nSpots = length(opts.subset);
end
opts.imageSize = job.ROI.cropSize;

% set up figure
figure(1); clf
fig_n=ceil(sqrt(nSpots));
fig_m=ceil(nSpots/fig_n);

% open movie
[~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata,0);
% read stack
img = kitReadImageStack(reader,job.metadata,opts.timePoint,opts.channel,job.ROI.crop,0);

% show each sister
for iSpot=opts.subset
    
    opts.spotID = iSpot;
    
    f = subplot(fig_m,fig_n,iSpot);
    
    noCoords = showSpot(job,img,opts);
    if noCoords
        cla(f);
    else
        plotTit = sprintf('Spot %i',iSpot);
        title(plotTit,'FontSize',10)
    end
end

warning(w);

end
   
%% SUB-FUNCTIONS
function noCoords = showSpot(job,img,opts)

% get coordinate system and plot channels
chan = opts.channel;
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,chan}(1:3);

% get sister information
tk = opts.spotID;
tp = opts.timePoint;

% get pixel resolution
pixelSize = job.metadata.pixelSize;
imageSize = size(img);

if isfield(job.dataStruct{chan},'tracks')
    % accumulate track information
    track = job.dataStruct{chan}.tracks(tk);
    startTime = track.seqOfEvents(1,1);
    endTime   = track.seqOfEvents(2,1);
    if tp < startTime || tp > endTime
        coords = nan(1,3);
    else
        coords = ...
            track.tracksCoordAmpCG(8*(tp-(startTime-1))-7:8*(tp-(startTime-1))-5);
        coords = coords + chrShift;
        coords = coords./pixelSize;
    end
else
    % accumulate initCoord information
    initCoord = job.dataStruct{chan}.initCoord(tp);
    coords = initCoord.allCoordPix(tk,1:3);
    chrShift = chrShift./pixelSize;
    coords = coords + chrShift;
end

% check whether any coordinates have been found, plot nothing if so
if any(isnan(coords))
    noCoords = 1;
    return
else
    noCoords = 0;
end

% calculate centre pixel
centrePxl = round(coords);

% max project over three z-slices around point
if opts.zProject == 1
    img = max(img, [], 3);
elseif opts.zProject == -1
    img = max(img(:,:,max(1,centrePxl(3)-2):min(centrePxl(3)+2,opts.imageSize(3))), [], 3);
else
    img = img(:,:,centrePxl(3));
end

% produce cropped image around track centre
xReg = [centrePxl(1)-ceil(0.5/pixelSize(1))+1 ...
           centrePxl(1)+ceil(0.5/pixelSize(1))+1];
yReg = [centrePxl(2)-ceil(0.5/pixelSize(2))+1 ...
           centrePxl(2)+ceil(0.5/pixelSize(2))+1];
xReg(1) = max(xReg(1),1); xReg(2) = min(xReg(2),imageSize(2));
yReg(1) = max(yReg(1),1); yReg(2) = min(yReg(2),imageSize(1));
imgCrpd = img(yReg(1):yReg(2),xReg(1):xReg(2));

% define contrast stretch and apply
irange=stretchlim(imgCrpd,opts.contrast);
imgCrpd = imadjust(imgCrpd, irange, []);

% correct coordinates to the cropped region
coords(:,1:2) = coords(:,1:2) - [xReg(1) yReg(1)] + 1;

% plot image
imshow(imgCrpd)

% if no coordinates to be plotted, skip the rest
if ~opts.plotCoords
    return
end

% plot coordinates
hold on
plot(coords(1),coords(2),...
    'Color','k','Marker','x','MarkerSize',10)
hold off

end
