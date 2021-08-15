function kitLabelSisters(job,varargin)
% KITLABELSISTERS Labels sisters in a given JOB.
%
%    KITLABELSISTERS(JOB,...) Labels sisters using lines over image of
%    movie from JOB for a single channel.
%
%    Options, defaults in {}:-
%
%    channel: 1, 2 or 3, where default is the {coordinate system channel}.
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
opts.timePoint = 1;
opts.transpose = 0;
opts.withinFig = 0;
% process user-defined options
opts = processOptions(opts, varargin{:});

%% Pre-processing

% convert coordinates to image into numbers
switch opts.coords
    case 'xy'
        opts.coords = [1 2];
        projectCoord = 3;
    case 'xz'
        opts.coords = [3 1];
        projectCoord = 1;
    case 'yz'
        opts.coords = [3 2];
        projectCoord = 2;
    otherwise
        error('Coordinates requested for imaging not recognised. Please provide either: ''xy'',''xz'' or ''yz''.');
end

% get the reader, and crop limits and size for the movie
if iscell(job)
    job = job{1};
    warning('Job provided is in cell format. Please ensure that you have provided a single job and not a full experiment.');
end
[md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'ROI',job.metadata,1);
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

% pre-define some structures
dS = job.dataStruct{opts.channel};
spotIDs = dS.sisterList(1).trackPairs(:,1:2);
hold on
% get sister coordinates
for iSis = 1:length(dS.sisterList)
    coords = NaN(2);
    for iSpot = 1:2
        iTrack = dS.tracks(spotIDs(iSis,iSpot));
        startT = iTrack.seqOfEvents(1,1);
        endT   = iTrack.seqOfEvents(2,1);
        if ismember(opts.timePoint,startT:endT)
            idx = dS.tracks(spotIDs(iSis,iSpot)).tracksFeatIndxCG(opts.timePoint-startT+1);
            if idx
                coords(iSpot,:) = dS.initCoord(opts.timePoint).allCoordPix(idx,opts.coords);
            end
        end
    end
    % plot
    line([coords(1,1) coords(2,1)],[coords(1,2) coords(2,2)],'Color','w');
end

end