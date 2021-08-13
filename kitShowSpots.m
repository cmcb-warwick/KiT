function imgCrpd = kitShowSpots(job,varargin)
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
%    cropHalfWidth: {0.5} number in um for size of box to crop spot to.
%    Will be overwritten if using marker 'i' to show radius of measurement.
%
%    marker: {'x'}, '+' or 'i'. Marker used to label spots. 'i' plots
%           circles with radius within which intensities were measured.
%
%    plotCoords: 0 or {1}. Whether or not to plot coordinates over images.
%
%    subset: {[]} or array of numbers. The spots to plot from a movie.
%           These numbers will refer to the track ID for each spot.
%
%    timePoint: {1} or number. Timepoint at which to plot spots.
%
%    withinFig: {0} or 1. Plot in within current figure environment or in 
%    a figure 1. 
%
%    zProject: {0}, 1 or -1. Whether or not to project in the z-direction.
%           -1 will project in the 5 z-slices surrounding the spot.
%
% Copyright (c) 2017 C. A. Smith

opts.channel = 1;
opts.contrast = [0.1 1];
opts.cropHalfWidth = 0.5;
opts.marker = 'x';
opts.noadjust = 0;
opts.plotCoords = 1;
opts.subset = [];
opts.timePoint = 1;
opts.useTrackID = 0;
opts.withinFig = 0;
opts.zProject = 0;
opts = processOptions(opts,varargin{:});

% suppress warnings
w = warning;
warning('off','all');

% check marker option
if ~ismember(opts.marker,'x+i')
    error('Marker option not recognised: %s',opts.marker);
end

% get important information
dataStruct = job.dataStruct{opts.channel};
if isempty(opts.subset)
%   if isfield(dataStruct,'tracks')
%     nSpots = length(dataStruct.tracks);
%     opts.subset = 1:nSpots;
%   else
    nSpots = size(dataStruct.initCoord(opts.timePoint).allCoord,1);
    opts.subset = 1:nSpots;
%   end
else
  nSpots = length(opts.subset);
end
if length(job.ROI)>1
    job.ROI = job.ROI(job.index);
    job = kitSaveJob(job);
end
opts.imageSize = job.ROI.cropSize;

% set up figure
if ~opts.withinFig
    figure(1); clf
end
fig_n=ceil(sqrt(nSpots));
fig_m=ceil(nSpots/fig_n);

% open movie
[~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata,0);
% read stack
img = kitReadImageStack(reader,job.metadata,opts.timePoint,opts.channel,job.ROI.crop,0);

% show each sister
subplot_ind = 0; %help to index over a subset of spots, rather than all
for iSpot=opts.subset
    subplot_ind=subplot_ind+1;
    opts.spotID = iSpot;
    
    f = subplot(fig_m,fig_n,subplot_ind);
    
    [noCoords, imgCrpd] = showSpot(job,img,opts);
    if noCoords
        set(f,'Visible','off');
    else
        plotTit = sprintf('Spot %i',iSpot);
        title(plotTit,'FontSize',10)
    end
end

warning(w);

end
   
%% SUB-FUNCTIONS
function [noCoords, imgCrpd] = showSpot(job,img,opts)

% get coordinate system and plot channels
chan = opts.channel;
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,chan}(1:3);

% set crop half-width
if strcmp(opts.marker,'i')
    maskRadius = job.options.intensity.maskRadius;
    cropHalfWidth = maskRadius*1.5;
else
    cropHalfWidth = opts.cropHalfWidth; %um
end

% get sister information
tk = opts.spotID;
tp = opts.timePoint;

% get pixel resolution
pixelSize = job.metadata.pixelSize;
imageSize = size(img);

% if isfield(job.dataStruct{chan},'tracks')
%     % accumulate track information
%     track = job.dataStruct{chan}.tracks(tk);
%     startTime = track.seqOfEvents(1,1);
%     endTime   = track.seqOfEvents(2,1);
%     if tp < startTime || tp > endTime
%         coords = nan(1,3);
%     else
%         coords = ...
%             track.tracksCoordAmpCG(8*(tp-(startTime-1))-7:8*(tp-(startTime-1))-5);
%         coords = coords + chrShift;
%         coords = coords./pixelSize;
%     end
% else
if ~opts.useTrackID
    % accumulate initCoord information
    initCoord = job.dataStruct{chan}.initCoord(tp);
    coords = initCoord.allCoordPix(tk,1:3);
else
        try
          coordsTMP = kitCoordsToImageCoords(job,chan,...
            job.dataStruct{chan}.trackList(tk).coords(tp,1:3),tp);
        catch
          coordsTMP = NaN(1,3);
        end
    coords = coordsTMP(1,:);
end
    chrShift = chrShift./pixelSize;
    coords = coords + chrShift;
%end

% check whether any coordinates have been found, plot nothing if so
if any(isnan(coords))
    noCoords = 1;
    imgCrpd = nan(2*ceil(cropHalfWidth/pixelSize(1)),2*ceil(cropHalfWidth/pixelSize(2)));
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
xReg = [centrePxl(1)-ceil(cropHalfWidth/pixelSize(1)) ...
           centrePxl(1)+ceil(cropHalfWidth/pixelSize(1))]; %previously was a +1 shift here. Seems this was not necessary?
yReg = [centrePxl(2)-ceil(cropHalfWidth/pixelSize(2)) ...
           centrePxl(2)+ceil(cropHalfWidth/pixelSize(2))];
xReg(1) = max(xReg(1),1); xReg(2) = min(xReg(2),imageSize(2));
yReg(1) = max(yReg(1),1); yReg(2) = min(yReg(2),imageSize(1));
imgCrpd = img(yReg(1):yReg(2),xReg(1):xReg(2));

% define contrast stretch and apply
irange=stretchlim(imgCrpd,opts.contrast);
if ~opts.noadjust
    imgCrpd = imadjust(imgCrpd, irange, []);
end

% correct coordinates to the cropped region
coords(:,1:2) = coords(:,1:2) - [xReg(1) yReg(1)] + 1;

% plot image
imshow(imgCrpd)

% if no coordinates to be plotted, skip the rest
if ~opts.plotCoords
    return
end

% plot coordinates.
hold on
if strcmp(opts.marker,'i')
    markerSize = ceil(job.options.intensity.maskRadius / pixelSize(1));
    markerSize = max(markerSize,0.5);
    drawCircle(coords(1),coords(2),markerSize,'w');
else
    markerSize = 10;
    plot(coords(1),coords(2),...
        'Color','k','Marker',opts.marker,'MarkerSize',markerSize)
end
hold off

end

function drawCircle(x,y,r,color)
% Draws circle.

% Estimate pixels in circumference.
c = 2*pi*r;
theta = linspace(0,2*pi,ceil(c));
cx = x + r*cos(theta);
cy = y + r*sin(theta);
plot(cx, cy, [color '-']);

end
