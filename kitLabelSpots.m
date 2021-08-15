function kitLabelSpots(job,varargin)
% KITLABELSPOTS Labels spots in a given JOB.
%
%    KITLABELSPOTS(JOB,...) Labels spot centres over image of movie from
%    JOB for a single channel.
%
%    Options, defaults in {}:-
%
%    channel: 1, 2 or 3, where default is the coordinate system channel.
%       Channel in which to overlay spots.
%
%    contrast: {[0.1 1]} or similar. Upper and lower contrast limits.
%       Values must be in range [0 1].
%       Tips: - Increase the lower limit to remove background noise.
%             - Decrease the upper limit to increase brightness.
%
%    coords: {'xy'}, 'xz' or 'yz'. Coordinate plane in which to show
%       images.
%
%    timePoint: {1} or positive integer. The time point of the movie at
%       which to show images.
%
%    transpose: {0} or 1. Whether or not to transpose images.
%
%    withinFig: {0} or 1. Whether or not to show images within the current
%       figure environment.
%
% Copyright (c) 2018 C. A. Smith

% define default options
opts.channel = job.options.coordSystemChannel;
opts.contrast = [0.1 1];
opts.coords = 'xy';
opts.marker = 'x';
opts.timePoint = 1;
opts.transpose = 0;
opts.withinFig = 0;
% process user-defined options
opts = processOptions(opts, varargin{:});

% check marker option
if ~ismember(opts.marker,'x+i')
    error('Marker option not recognised: %s',opts.marker);
end

%% Pre-processing

% convert coordinates to image into numbers
switch opts.coords
    case 'xy'
        opts.coords = [1 2];
        projectCoord = 3;
    case 'xz'
        opts.coords = [3 1];
        projectCoord = 1; %x and y switched in coords
    case 'yz'
        opts.coords = [3 2];
        projectCoord = 2; %x and y switched in coords
    otherwise
        error('Coordinates requested for imaging not recognised. Please provide either: ''xy'',''xz'' or ''yz''.');
end

% get the reader, and crop for the movie
if iscell(job)
    job = job{1};
    warning('Job provided is in cell format. Please ensure that you have provided a single job and not a full experiment.');
end
[md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
crop = job.ROI.crop;

%% Process the image

% get image stack for this timepoint
img = kitReadImageStack(reader, md, opts.timePoint, opts.channel, crop, 0);

% project the image in the projection coordinate
img = max(img,[],projectCoord);
img = squeeze(img);
if opts.transpose
    img = img';
end

% contrasting
irange = stretchlim(img,opts.contrast);
img = imadjust(img, irange, []);

% plotting

if ~opts.withinFig
    figure(1)
    clf
end
imshow(img);

%% Processing coordinate information

% get spot coordinates
dS = job.dataStruct{opts.channel};
coords = dS.initCoord(opts.timePoint).allCoordPix(:,opts.coords);
hold on
if strcmp(opts.marker,'i')
    markerSize = ceil(job.options.intensity.maskRadius / job.metadata.pixelSize(1));
    for i = 1:size(coords,1)
        drawCircle(coords(i,1),coords(i,2),markerSize,'w');
    end
else
    markerSize = 10;
    scatter(coords(:,1),coords(:,2),['w' opts.marker])
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