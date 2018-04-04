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
%
% Copyright (c) 2017 C. A. Smith

opts.channel = 1;
opts = processOptions(opts,varargin{:});

% some variables
gridsep = 2;
imgHalfWidth = 0.5; %um
opts.bgcol = dot([0.94 0.94 0.94],[0.2989 0.5870 0.1140]);

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
    nROIs = job.nROIs;
elseif length(job.ROI)>1
    nROIs = length(job.ROI);
else
    nROIs = NaN;
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
if isfield(dS,'rawData')
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
figtit = sprintf('Spot filtering: Image %i of %i, channel %i',job.index,nROIs,opts.channel);
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
   
%% SUB-FUNCTIONS
function [imgCrpd,coords] = indivSpot(img,coords,opts)

% get pixel resolution
imageSize = size(img);
imgw = opts.imgHalfWidth;
bgcol = opts.bgcol;

% predefine cropped image
imgCrpd = ones(2*imgw+1)*bgcol;

% calculate centre pixel
centrePxl = round(coords);

% max project over three z-slices around point
img = max(img(:,:,max(1,centrePxl(3)-2):min(centrePxl(3)+2,opts.imageSize(3))), [], 3);

% produce cropped image around track centre
xReg = [max(centrePxl(1)-imgw+1,1) min(centrePxl(1)+imgw+1,imageSize(2))];
yReg = [max(centrePxl(2)-imgw+1,1) min(centrePxl(2)+imgw+1,imageSize(1))];
imgCrpd(1:diff(yReg)+1,1:diff(xReg)+1) = img(yReg(1):yReg(2),xReg(1):xReg(2));

% define contrast stretch and apply
irange=stretchlim(imgCrpd,[0.1 1]);
imgCrpd = imadjust(imgCrpd, irange, []);

% correct coordinates to the cropped region
coords(:,1:2) = coords(:,1:2) - [xReg(1) yReg(1)];

end