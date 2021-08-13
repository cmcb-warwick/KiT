function imgDims = griddedSpots(job,varargin)
% IMGDIMS = GRIDDEDSPOTS Plots images of each spot in a grid for use by
% kitFilterSpots.
%
%    IMGDIMS = GRIDDEDSPOTS(JOB,...) Plots coordinates over images of each
%    spot localised in a given channel, outputting dimensions and positions
%    of each element of the grid.
%
%    Options, defaults in {}:-
%
%    channel: {1}, 2, 3 or 4. Which channel to show.
%
%    rawData: 0 or {1}. Whether or not to show raw data.
%
%
% Copyright (c) 2017 C. A. Smith

opts.channel = 1; %default options
opts.rawData = 1;
opts = processOptions(opts,varargin{:});

% some variables
gridsep = 2;
imgHalfWidth = 0.5; %um
opts.bgcol = dot([0.94 0.94 0.94],[0.2989 0.5870 0.1140]); %background colour

% calculate imgWidth in pixels
pixelSize = job.metadata.pixelSize;
opts.imgHalfWidth = ceil(imgHalfWidth/pixelSize(1));
imgWidth = 2*opts.imgHalfWidth+1;

% suppress warnings
w = warning;
warning('off','all');

% get important information
dS = job.dataStruct{opts.channel};
nSpots = size(dS.initCoord.allCoord,1);
opts.imageSize = job.ROI.cropSize;
if isfield(job,'nROIs')
    nROIs = [' of ' num2str(job.nROIs)];
else
    nROIs = '';
end

% set up figure
figure(1); clf
fig_n=ceil(sqrt(nSpots));
fig_m=ceil(nSpots/fig_n);

% open movie and read stack
if length(job.ROI)>1
    [~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI(job.index).movie),'ROI');
    img = kitReadImageStack(reader,job.metadata,1,opts.channel,job.ROI(job.index).crop,0);
else
    [~,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'ROI');
    img = kitReadImageStack(reader,job.metadata,1,opts.channel,job.ROI.crop,0);
end

% make empty gridded image
gridw = fig_m*imgWidth + (fig_m+1)*gridsep;
gridh = fig_n*imgWidth + (fig_n+1)*gridsep;
gridImg = opts.bgcol*ones(gridh,gridw);

% get all coordinate
if isfield(dS,'rawData') && opts.rawData
  iC = dS.rawData.initCoord;
else
  iC = dS.initCoord;
end
allCoords = iC.allCoordPix(:,1:3);

% form all spot positions
spotpos = 1:nSpots;
spotpos = spotpos(:);
spotpos = [mod(spotpos,fig_m) ceil(spotpos./fig_m)];
spotpos(spotpos(:,1)==0,1) = fig_m;

% preset allRnge array
rnge = [gridsep+(gridsep+imgWidth).*(spotpos(:,1)-1)+1 ...
        gridsep+(gridsep+imgWidth).*(spotpos(:,2)-1)+1];

% show each sister
for iSpot=1:nSpots
    
    % check whether any coordinates have been found, do nothing if so
    coords = allCoords(iSpot,:);
    if any(isnan(coords))
        allCoords(iSpot,:) = NaN;
        continue
    end
    
    % calculate the pixels in which to push new image
    [gridImg(rnge(iSpot,2):rnge(iSpot,2)+imgWidth-1, ...
        rnge(iSpot,1):rnge(iSpot,1)+imgWidth-1),coords] ...
        = indivSpot(img,coords,opts);
    
    % calculate position of coordinates to be plotted
    allCoords(iSpot,1:2) = coords(:,1:2)+rnge(iSpot,:);
    
end

% plot the full image and coordinates
imshow(gridImg,'Border','tight');
hold on
scatter(allCoords(:,1),allCoords(:,2),15*fig_m,'k','x')
figtit = sprintf('Spot filtering: Image %i%s, channel %i',job.index,nROIs,opts.channel);
set(gcf,'Resize','off','Name',figtit,'Units','characters',...
    'Position',[70 35 80 50],'NumberTitle','off');
movegui(gcf,'center');

% save image dimensions and positions
imgDims = [rnge repmat(imgWidth,nSpots,1)];

% close the reader
close(reader);

% reset warnings
warning(w);

end
   
