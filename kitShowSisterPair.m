function kitShowSisterPair(job,varargin)
% KITSHOWSISTERPAIR Plots image of movie with coordinates for a given
% sister pair.
%
%    KITSHOWSISTERPAIR(JOB,...) Plots coordinates over the image of a given
%    sister pair in a given channel at a given timepoint.
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
%    newFig: {0} or 1. Whether or not to show the sister pair in a new
%           figure.
%
%    sisterPair: {1} or number. Sister pair within JOB being plotted.
%
%    timePoint: {1} or number. Timepoint at which to plot the sister pair.
%
%    title: {[]} or a string. The title for the image. When left empty, the
%           title will reflect the sister number and time point chosen.
%
%    transpose: {0} or 1. Whether to transpose the image.
%
%    withinFig: {0} or 1. Whether or not to show the sister pair within the
%           current figure environment.
%
%    zoomScale: {1} or number. Magnification of a zoomed image as a
%           proportion of the default, so that 0.5 will zoom out, and 1.5
%           will zoom in.
%
%    zoom: 0 or {1}. Whether or not to zoom into the specific sister
%           pair. A value of 0 plots the whole cell, with just the
%           coordinates of the sister pair plotted.
%
%    zProject: 0, 1 or {-1}. Whether or not to project in the z-direction.
%           -1 will project in the 5 z-slices surrounding the sister pair.
%
% Copyright (c) 2017 C. A. Smith

if nargin<1
  error('Must supply a job.');
end

% set default options
opts.contrast = [0.1 1];
opts.channel = 1;
opts.newFig = 0;
opts.sisterPair = 1;
opts.timePoint = 1;
opts.title = [];
opts.transpose = 0;
opts.withinFig = 0;
opts.zoomScale = 1;
opts.zoom = 1;
opts.zProject = -1;

% process options
opts = processOptions(opts, varargin{:});

%% Image and coordinate acquisition

% open movie
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata,0);

% get coordinate system and plot channels
chan = opts.channel;
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,chan}(1:3);

% get sister information
sisPair = opts.sisterPair;
sisterList = job.dataStruct{chan}.sisterList;

if sisPair > length(sisterList)
    kitLog('Sister pair ID provided too large. Only have %i sisters in movie. Quitting.',length(sisterList))
    return
end

% get track information
trackIDs = sisterList(1).trackPairs(sisPair,1:2);
timePoint = opts.timePoint;

% get pixel resolution
pixelSize = job.metadata.pixelSize;

coords = nan(2,3); %coords x sister

% accumulate track information by channel and sister
for iSis = 1:2
    tk = trackIDs(iSis);
    track = job.dataStruct{chan}.tracks(tk);

    startTime = track.seqOfEvents(1,1);
    endTime   = track.seqOfEvents(2,1);
    if timePoint < startTime || timePoint > endTime
        coords(iSis,:) = nan(1,3);
    else
        coords(iSis,:) = ...
            track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
        coords(iSis,:) = coords(iSis,:) + chrShift;
        coords(iSis,:) = coords(iSis,:)./pixelSize;
    end
end

if sum(isnan(coords(:))) == 6
    kitLog('No coordinates found for sister pair %i at time point %i. Ignoring.',sisPair,timePoint);
    return
end

% calculate pair centre and convert to pixels
centrePxl = nanmean(coords);
centrePxl = round(centrePxl);


%% Image production

% read stack
img = kitReadImageStack(reader,md,timePoint,chan,job.ROI.crop,0);
% max project over three z-slices around point
if opts.zProject == 1
    img = max(img, [], 3);
elseif opts.zProject == -1
    img = max(img(:,:,centrePxl(3)-2:centrePxl(3)+2), [], 3);
else
    img = img(:,:,centrePxl(3));
end
if opts.transpose
    img = img';
end

% produce cropped image around track centre
if opts.zoom
    
    xReg = [centrePxl(1)-ceil(opts.zoomScale*(2/pixelSize(1)))+1 ...
               centrePxl(1)+ceil(opts.zoomScale*(2/pixelSize(1)))+1];
    yReg = [centrePxl(2)-ceil(opts.zoomScale*(2/pixelSize(2)))+1 ...
               centrePxl(2)+ceil(opts.zoomScale*(2/pixelSize(2)))+1];
           
    if opts.transpose
        imgCrpd = img(xReg(1):xReg(2),yReg(1):yReg(2));
    else
        imgCrpd = img(yReg(1):yReg(2),xReg(1):xReg(2));
    end
    
    % define contrast stretch of cropped region and apply
    irange=stretchlim(imgCrpd,opts.contrast);
    imgCrpd = imadjust(imgCrpd, irange, []);
    
    % correct coordinates to the cropped region
    coords(:,1:2) = coords(:,1:2) - repmat([xReg(1) yReg(1)],2,1) + 1;

else
    % define contrast stretch of full image and apply
    irange=stretchlim(img,opts.contrast);
    img = imadjust(img, irange, []);
end


%% Producing the figure

% prepare figure environment
if opts.newFig
    figure;
elseif ~opts.withinFig
    figure(1)
    clf
end

% plot image
if opts.zoom
    imshow(imgCrpd)
else
    imshow(img)
end

% plot coordinates
hold on
for iSis = 1:2
    if opts.transpose
        plot(coords(iSis,2),coords(iSis,1),...
            'Color','k','Marker','x','MarkerSize',15)
    else
        plot(coords(iSis,1),coords(iSis,2),...
            'Color','k','Marker','x','MarkerSize',15)
    end
end

% aesthetics
if isempty(opts.title)
  opts.title = sprintf('Sister pair %i, time point %i',opts.sisterPair,opts.timePoint);
end
title(opts.title,'FontSize',20)

hold off

end