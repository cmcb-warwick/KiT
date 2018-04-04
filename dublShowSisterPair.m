function dublShowSisterPair(job,varargin)
% DUBLSHOWSISTERPAIR Plots image of dual-channel movie with coordinates
% for a given sister pair.
%
%    DUBLSHOWSISTERPAIR(JOB,...) Plots coordinates in two channels over
%    the movie image for a given sister pair at a given timepoint.
%
%    Options, defaults in {}:-
%
%    channelMap: {[2 1 3]} or some perturbation. Order in which the titled
%    channels are presented in the figure, where typically:
%    1=red, 2=green, 3=blue.
%
%    contrast: {[0.1 1]} or similar structure. A 1x3 cell of pairs of
%           numbers for use in image contrast.
%
%    newFig: {0} or 1. Whether or not to show the sister pair in a new
%           figure.
%
%    plotChannels: {[1 2]} or some subset of [1 2 3]. Vector of channels
%           for plotting, where typically:
%               1=red, 2=green, 3=blue.
%
%    plotTitles: {'mNeonGreen','tagRFP','JF-646'} or similar. Titles of
%           each plot denoting which channel is which.
%
%    sisterPair: {1} or number. Sister pair within JOB being plotted.
%
%    subpixelate: {9} or odd number. Number of pixels of accuracy within
%           which to correct for chromatic shift.
%
%    timePoint: {1} or number. Timepoint at which to plot the sister pair.
%
%    transpose: {0} or 1. Whether to tranpose the image.
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
% Copyright (c) 2014 C. A. Smith

if nargin<1
  error('Must supply a job.');
end

% set default options
opts.channelMap = [2 1 3]; % green, red, blue
opts.contrast = repmat({[0.1 1]},1,3);
opts.newFig = 0;
opts.plotChannels = 1:2;
opts.plotTitles = {'mNeonGreen','tagRFP','JF-646'};
opts.sisterPair = 1;
opts.subpixelate = 9;
opts.timePoint = 1;
opts.transpose = 0;
opts.zoomScale = 1;
opts.zoom = 1;
opts.zProject = -1;

% process options
opts = processOptions(opts, varargin{:});
while length(opts.plotTitles) < 3
    opts.plotTitles = [opts.plotTitles,''];
end

%% Image and coordinate acquisition

% open movie
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);

% get coordinate system and plot channels
coordSysChan = job.options.coordSystemChannel;
plotChans = opts.plotChannels;
nChans = length(plotChans);

% get sister information
sisPair = opts.sisterPair;
sisterList = job.dataStruct{coordSysChan}.sisterList;

% get track information
trackIDs = sisterList(1).trackPairs(sisPair,1:2);
timePoint = opts.timePoint;

% if only one channel, run the single channel version
if nChans == 1
    kitLog('Only one channel being shown. Running kitShowSisterPair instead.');
    kitShowSisterPair(job,'channel',plotChans,'contrast',opts.contrast,'newFig',opts.newFig,...
        'sisterPair',sisPair,'timePoint',timePoint,'title',opts.plotTitles{1},'transpose',opts.transpose,...
        'withinFig',0,'zoomScale',opts.zoomScale,'zoom',opts.zoom,'zProject',opts.zProject);
    return
end

% get pixel resolution
pixelSize = job.metadata.pixelSize;

% accumulate track information by channel and sister
trackCoord = nan(nChans,3,2);
for c = plotChans
    for iSis = 1:2
        tk = trackIDs(iSis);
        track = job.dataStruct{c}.tracks(tk);
        
        startTime = track.seqOfEvents(1,1);
        if timePoint < startTime
            trackCoord(c,:,iSis) = nan(1,3);
        else
            trackCoord(c,:,iSis) = ...
                track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
        end
    end
end

% calculate pair centre - check whether coordinate system is reference
if ismember(coordSysChan,plotChans)
    centrePoint = nanmean(trackCoord(coordSysChan,:,:),3);
else
    % if not, then use first available channel
    coordSysChan = min(plotChans);
    centrePoint = nanmean(trackCoord(coordSysChan,:,:),3);
end
% convert to pixels
centrePoint = centrePoint./pixelSize;
centrePxl = round(centrePoint);

mapChans = opts.channelMap;
if opts.transpose
    dims = [2 1];
else
    dims = [1 2];
end

%% RGB image production

% predesignation of images
rgbImg = zeros([job.ROI.cropSize(dims), 3]);
rgbImgCS = zeros([job.ROI.cropSize(dims)*opts.subpixelate, 3]);
            
if opts.zoom
    % calculate spread about the centre pixels
    cropSpread = opts.zoomScale*(2./pixelSize);
	cropSpread = ceil(cropSpread);
    rgbCrpd= zeros([2*cropSpread(1)+1  2*cropSpread(2)+1  3]);
    rgbCrpdCS= zeros([opts.subpixelate*(2*cropSpread(1)+1) ...
                opts.subpixelate*(2*cropSpread(2)+1) ...
                3]);
end

% produce raw image
for c = plotChans
    
    % read stack
    img = kitReadImageStack(reader,md,timePoint,c,job.ROI.crop,0);
    
    % max project over 5 z-slices around point
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
    rgbImg(:,:,mapChans(c)) = img;
    
end

% produce cropped image around pair centre
if opts.zoom
    
    % non-chromatic shifted regions
    xReg = [centrePxl(1)-cropSpread(1) centrePxl(1)+cropSpread(1)];
    yReg = [centrePxl(2)-cropSpread(2) centrePxl(2)+cropSpread(2)];
    
    for c = plotChans
        if opts.transpose
            rgbCrpd(:,:,mapChans(c)) = rgbImg(xReg(1):xReg(2),yReg(1):yReg(2),mapChans(c));
        else
            rgbCrpd(:,:,mapChans(c)) = rgbImg(yReg(1):yReg(2),xReg(1):xReg(2),mapChans(c));
        end
    end
    
    % produce chromatic shifted image
    first = 1;
    for c = setdiff(plotChans,coordSysChan) 
      if first

        [img1,img2] = ... 
            chrsComputeCorrectedImage(rgbCrpd(:,:,mapChans(coordSysChan)),rgbCrpd(:,:,mapChans(c)),job.options.chrShift.result{coordSysChan,c}, ...
            'subpixelate',opts.subpixelate);

        % give the subpixelated images to the image structure
        rgbCrpdCS(:,:,mapChans(coordSysChan)) = img1;
        rgbCrpdCS(:,:,mapChans(c)) = img2;

        first = 0;

      else
        [~,img2] = ...
            chrsComputeCorrectedImage(rgbCrpd(:,:,mapChans(coordSysChan)),rgbCrpd(:,:,mapChans(c)),job.options.chrShift.result{coordSysChan,c}, ...
            'subpixelate',opts.subpixelate);

        % give the new subpixelated images to the image structure
        rgbCrpdCS(:,:,mapChans(c)) = img2;

      end
    end
    
    % define contrast stretch for shifted cropped, and apply
    for c = plotChans
        irange = stretchlim(rgbCrpdCS(:,:,mapChans(c)),opts.contrast{c});
        rgbCrpdCS(:,:,mapChans(c)) = imadjust(rgbCrpdCS(:,:,mapChans(c)), irange, []);
    end
    
else
    
    % produce chromatic shifted image
    first = 1;
    for c = setdiff(plotChans,coordSysChan) 
      if first
        
        [img1,img2] = ... 
            chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(c)),job.options.chrShift.result{coordSysChan,c}, ...
            'subpixelate',opts.subpixelate);
        
        % give the subpixelated images to the image structure
        rgbImgCS(:,:,mapChans(coordSysChan)) = img1;
        rgbImgCS(:,:,mapChans(c)) = img2;
        
        first = 0;
    
      else
        [~,img2] = ...
            chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(c)),job.options.chrShift.result{coordSysChan,c}, ...
            'subpixelate',opts.subpixelate);
        
        % give the new subpixelated images to the image structure
        rgbImgCS(:,:,mapChans(c)) = img2;
        
      end
    end
    
    % define contrast stretch for cropped images, and apply
    for c = plotChans
        irange = stretchlim(rgbImgCS(:,:,mapChans(c)),opts.contrast{c});
        rgbImgCS(:,:,mapChans(c)) = imadjust(rgbImgCS(:,:,mapChans(c)), irange, []);
    end
    
end


%% Find coordinates to plot for each RGB image
    
% transform track coordinates to centre of cropped region
coord = trackCoord(plotChans,:,:);

% adjust to pixels
for i = 1:2
    coord(:,:,i) = coord(:,:,i)./repmat(pixelSize,nChans,1);
end
coordCS = coord*opts.subpixelate;    

% correct coordinates for region position
if opts.zoom
    coord(:,1,:) = coord(:,1,:) - (xReg(1)+1);
    coord(:,2,:) = coord(:,2,:) - (yReg(1)-1);
    coordCS(:,1,:) = coordCS(:,1,:) - xReg(1)*opts.subpixelate + (opts.subpixelate+1)/2; % for subpix=9, +5; subpix=3, +2; subpix=5, +3 
    coordCS(:,2,:) = coordCS(:,2,:) - yReg(1)*opts.subpixelate + (opts.subpixelate+1)/2;
end

%% Producing the figure

if nChans == 1
    plotImg = rgb2gray(rgbImg);
    if opts.zoom; plotImgCrpd = rgb2gray(rgbCrpd); end
    plotCoord = coord;
else
    plotImg = rgbImgCS;
    if opts.zoom; plotImgCrpd = rgbCrpdCS; end
    plotCoord = coordCS;
end 

% prepare figure environment
if opts.newFig
    figure
else
    figure(1)
end
clf

C = [ 0  1  0;
      1  0  0;
      0  0  1];
plotStyle = ['x' '+' 'o'];

if opts.zoom
    
    % plot individual channels
    bigImgInd = [];
    if nChans > 1
        for c = 1:nChans
            subplot(nChans,nChans+1,c*(nChans+1))
            hold on
            imshow(plotImgCrpd(:,:,mapChans(c)))
            title(opts.plotTitles{c},'FontSize',20,'Color',C(plotChans(c),:))
            bigImgInd = [bigImgInd, (c-1)*(nChans+1)+1 : c*(nChans+1)-1 ];
        end
    end
    % plot dual image
    if isempty(bigImgInd)
        bigImgInd = 1;
    end
    subplot(nChans,nChans+1,bigImgInd)
    imshow(plotImgCrpd)
    title('All channels','FontSize',20)

    hold on
    % plot coordinates
    for c = 1:nChans
        for i = 1:2
            subplot(nChans,nChans+1,bigImgInd)
            if opts.transpose
                plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                    'Color',C(plotChans(c),:),'Marker',plotStyle(plotChans(c)),'MarkerSize',15)
            else
                plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                    'Color',C(plotChans(c),:),'Marker',plotStyle(plotChans(c)),'MarkerSize',15)
            end
            if nChans > 1
                subplot(nChans,nChans+1,c*(nChans+1))
                if opts.transpose
                    plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                        'Color','k','Marker',plotStyle(plotChans(c)),'MarkerSize',15)
                else
                    plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                        'Color','k','Marker',plotStyle(plotChans(c)),'MarkerSize',15)
                end
            end
        end
    end
    
else
    % plot image
    imshow(plotImg)
    
    hold on
    % plot tracked channel's coordinates (points too close together on full image)
    for c = 1:nChans
        for i = 1:2
            if opts.transpose
                plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                    'Color',C(plotChans(c),:),'Marker',plotStyle(plotChans(c)),'MarkerSize',15)
            else
                plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                    'Color',C(plotChans(c),:),'Marker',plotStyle(plotChans(c)),'MarkerSize',15)
            end
        end
    end
end

hold off

end