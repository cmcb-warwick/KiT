function kitShowSisters(job,varargin)
% KITSHOWSISTERS Plots images of each sister pair in a movie.
%
%    KITSHOWSISTERS(JOB,...) Plots coordinates over images of each
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
%    plotCoords: 0 or {1}. Whether or not to plot coordinates over images.
%
%    subset: {[]} or array of numbers. The sisters to plot from a movie.
%
%    timePoint: {1} or number. Timepoint at which to plot sister pairs.
%
%    zProject: 0, 1 or {-1}. Whether or not to project in the z-direction.
%           -1 will project in the 5 z-slices surrounding the sister pair.
%
% Copyright (c) 2017 C. A. Smith

opts.channel = 1;
opts.contrast = [0.1 1];
opts.plotCoords = 1;
opts.subset = [];
opts.timePoint = 1;
opts.zProject = -1;
opts = processOptions(opts,varargin{:});

% get important information
dataStruct = job.dataStruct{opts.channel};
if isempty(opts.subset)
  nSisters = length(dataStruct.sisterList);
  opts.subset = 1:nSisters;
else
  nSisters = length(opts.subset);
end
opts.imageSize = job.ROI.cropSize;

% set up figure
figure(1); clf
fig_n=ceil(sqrt(nSisters));
fig_m=ceil(nSisters/fig_n);

% open movie
[~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata,0);
% read stack
img = kitReadImageStack(reader,job.metadata,opts.timePoint,opts.channel,job.ROI.crop,0);

% show each sister
for iSis=opts.subset
    
    opts.sisterPair = iSis;
    
    f = subplot(fig_m,fig_n,iSis);
    
    noCoords = showSisterPair(job,img,opts);
    if noCoords
        cla(f);
    else
        plotTit = sprintf('Sis %i',iSis);
        title(plotTit,'FontSize',10)
    end
end

end
   
%% SUB-FUNCTIONS
function noCoords = showSisterPair(job,img,opts)

% get coordinate system and plot channels
chan = opts.channel;
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,chan}(1:3);

% get sister information
sisID = opts.sisterPair;
sisterList = job.dataStruct{chan}.sisterList;

% get track information
trackIDs = sisterList(1).trackPairs(sisID,1:2);
tP = opts.timePoint;

% get pixel resolution
pixelSize = job.metadata.pixelSize;

coords = nan(2,3); %coords x sister

% accumulate track information by channel and sister
for iSis = 1:2
    tk = trackIDs(iSis);
    track = job.dataStruct{chan}.tracks(tk);

    startTime = track.seqOfEvents(1,1);
    endTime   = track.seqOfEvents(2,1);
    if tP < startTime || tP > endTime
        coords(iSis,:) = nan(1,3);
    else
        coords(iSis,:) = ...
            track.tracksCoordAmpCG(8*(tP-(startTime-1))-7:8*(tP-(startTime-1))-5);
        coords(iSis,:) = coords(iSis,:) + chrShift;
        coords(iSis,:) = coords(iSis,:)./pixelSize;
    end
end

% check whether any coordinates have been found, plot nothing if so
if sum(isnan(coords(:))) == 6
    noCoords = 1;
    return
else
    noCoords = 0;
end

% calculate pair centre and convert to pixels
centrePxl = nanmean(coords);
centrePxl = round(centrePxl);

% max project over three z-slices around point
if opts.zProject == 1
    img = max(img, [], 3);
elseif opts.zProject == -1
    img = max(img(:,:,max(1,centrePxl(3)-2):min(centrePxl(3)+2,opts.imageSize(3))), [], 3);
else
    img = img(:,:,centrePxl(3));
end

% produce cropped image around track centre
xReg = [centrePxl(1)-ceil(1/pixelSize(1))+1 ...
           centrePxl(1)+ceil(1/pixelSize(1))+1];
yReg = [centrePxl(2)-ceil(1/pixelSize(2))+1 ...
           centrePxl(2)+ceil(1/pixelSize(2))+1];
imgCrpd = img(yReg(1):yReg(2),xReg(1):xReg(2));

% define contrast stretch and apply
irange=stretchlim(imgCrpd,opts.contrast);
imgCrpd = imadjust(imgCrpd, irange, []);

% correct coordinates to the cropped region
coords(:,1:2) = coords(:,1:2) - repmat([xReg(1) yReg(1)],2,1) + 1;

% plot image
imshow(imgCrpd)

% if no coordinates to be plotted, skip the rest
if ~opts.plotCoords
    return
end

% plot coordinates
hold on
for iSis = 1:2
    plot(coords(iSis,1),coords(iSis,2),...
        'Color','k','Marker','x','MarkerSize',15)
end
hold off

end
