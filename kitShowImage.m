function kitShowImage(job,varargin)
% KITSHOWIMAGE Plots image of a multi-channel movie.
%
%    KITSHOWIMAGE(JOB,...) Shows cropped image of movie from JOB.
%
%    Options, defaults in {}:-
%
%    chanOrder: {'rgb'} or other permutation of r (red), g (green) and b
%       (blue). Order in which channels provided in the job file are
%       represented as colours in the rgb image.
%
%    contrast: {{[0.1 1],[0.1 1],[0.1 1]}}, 'help', or similar. Upper and
%       lower contrast limits for each channel. Values must be in range
%       [0 1]. 'help' outputs the minimum and maximum intensity values as 
%       guidance, then requests values from the user.
%       Tips: - Increase the lower limit to remove background noise.
%             - Decrease the upper limit to increase brightness.
%
%    coords: {'xy'}, 'xz' or 'yz'. Coordinate plane in which to show
%       images.
%
%    crop: 0, {1} or -1. Whether or not to show cropped images, as defined
%       by ROI.crop. -1 shows a full image with the cropped region
%       annotated.
%
%    jobsetMovie: {[]} or positive integer. The number movie to use when
%       providing a jobset rather than job structure.
%
%    imageChans: {[1 2]} or subset of [1,2,3]. Channels in which to show
%       images.
%
%    projectionRange: {[]}, 'help', or subset of other coordinate's full
%       range. The range of pixels over which to project the third
%       coordinate not given by 'coords'. 'help' outputs the possible
%       values the range can take, then requests a range from the user.
%
%    scaleBarLabel: {0} or 1. Whether or not to label the scale bar's
%       length.
%
%    scaleBarSize: {3} or distance in µm. Length of the scale bar. Scale
%       bars are overlaid at the south west (bottom right) of the image.
%
%    subpixelate: {9} or positive integer. The number by which to divide
%       pixels in order to allow accurate chromatic shift of images.
%
%    textNorthWest, textNorthEast and textSouthEast: {[]} or a string.
%       String labels to be overlaid on the image in the designated corner
%       of the image. To separate text into multiple lines in a given
%       corner, provide each line as a string in a different element of a
%       cell, e.g. {'mNeonGreen-','CENP-A'}.
%
%    timePoint: {1} or positive integer. The time point of the movie at
%       which to show images.
%
%    transpose: {0} or 1. Whether or not to transpose images.
%
%    withinFig: {0} or 1. Whether or not to show images within the current
%       figure environment.
%
% Copyright (c) 2017 C. A. Smith

% define default options
opts.chanOrder = 'grb';
opts.contrast = repmat({[0.1 1]},1,3);
opts.coords = 'xy';
opts.crop = 1; % use -1 to show crop on a full image
opts.jobsetMovie = [];
opts.imageChans = [1 2];
opts.projectionRange = [];
opts.scaleBarLabel = 0;
opts.scaleBarSize = 3; % in µm, length of scale bar
opts.subpixelate = 9;
opts.textNorthWest = [];
opts.textNorthEast = [];
opts.textSouthEast = [];
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
    case 'xz'
        opts.coords = [1 3];
    case 'yz'
        opts.coords = [2 3];
    otherwise
        error('Coordinates requested for imaging not recognised. Please provide either: ''xy'',''xz'' or ''yz''.');
end

% convert channel order into numbers
if length(opts.chanOrder) < length(opts.imageChans)
    error('Number of channels in channel order (n=%i) does not reflect the number of channels being imaged (n=%i).',length(opts.chanOrder),length(opts.imageChans));
else
    chanOrder = [3 3 3];
    colours = {'r','g','b'};
    for iChan = 3:-1:1
        tempLoc = find(opts.chanOrder == colours{iChan});
        if ~isempty(tempLoc)
            chanOrder(tempLoc) = iChan;
        end
    end
    if sum(chanOrder == 3) > 1
        error('More channels required in channel order. Please provide at least 2.');
    end
end

% get the reader, and crop limits and size for the movie
if isempty(opts.jobsetMovie)
    if iscell(job)
        job = job{1};
        warning('Job provided is in cell format. Please ensure that you have provided a single job and not a full experiment.');
    end
    jobID = job.index;
    [md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);
    if opts.crop==1
      crop = job.ROI.crop;
      cropSize = job.ROI.cropSize;
    else
      crop = [];
      cropSize = md.frameSize;
    end
else
    jobID = opts.jobsetMovie;
    if isfield(job,'metadata')
      [md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata{jobID});
    else
      [md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
    end
    if opts.crop==1
      crop = job.ROI(jobID).crop;
      cropSize = job.ROI(jobID).cropSize;
    else
      crop = [];
      cropSize = md.frameSize;
    end
end

% get required metadata
coordSysChan = job.options.coordSystemChannel;
pixelSize = md.pixelSize(1:3);
chrShift = job.options.chrShift.result;

% subpixelate value
if length(opts.imageChans)==1
    opts.subpixelate = 1;
end

% provide guidance on projection range if requested
projectCoord = setdiff(1:3,opts.coords);
if strcmp(opts.projectionRange,'help')
    % get maximum crop size for third coordinate
    maxOtherCoordRange = cropSize(projectCoord);
    % output possible range 
    fprintf('\nPROJECTION RANGE GUIDANCE:\nRange in third coordinate: [1 %i].\n',maxOtherCoordRange);
    % request input from user
    userRange = input('Please provide range over which to project third coordinate: ');
    while userRange(1)<1 || userRange(2)>maxOtherCoordRange || userRange(1)>userRange(2) 
        userRange = input('These values fall outside the permitted range. Please try again: ');
    end
    opts.projectionRange = userRange(1):userRange(2);
% otherwise check whether one has been provided
elseif isempty(opts.projectionRange)
    opts.projectionRange = 1:cropSize(projectCoord);
end

%% Create image structures

% produce structure to hold images, dependent on transposing
if opts.transpose
    rgbImg = zeros([cropSize(fliplr(opts.coords)), 3]);
    rgbImgShift = zeros([cropSize(fliplr(opts.coords))*opts.subpixelate, 3]);
else
    rgbImg = zeros([cropSize(opts.coords), 3]);
    rgbImgShift = zeros([cropSize(opts.coords)*opts.subpixelate, 3]);
end

% get images for each channel
if strcmp(opts.contrast,'help')
    fprintf('\nCONTRAST GUIDANCE:\n');
end
for iChan = opts.imageChans
    
    % get image stack for this timepoint
    img = kitReadImageStack(reader, md, opts.timePoint, iChan, crop, 0);
    
    % project the image in the third coordinate
    switch projectCoord
        case 1
            img = img(opts.projectionRange,:,:);
        case 2
            img = img(:,opts.projectionRange,:);
        case 3
            img = img(:,:,opts.projectionRange);
    end
    img = max(img,[],projectCoord);
    % reshape the image to ensure is in the correct form
    img = reshape(img, size(img,opts.coords(1)), size(img,opts.coords(2)));
    if opts.transpose
        img = img';
    end
    
    % contrasting
    if iscell(opts.contrast)
        
        irange(iChan,:) = stretchlim(img,opts.contrast{iChan});
        
    elseif strcmp(opts.contrast,'help')
        
        % get maximum and minimum intensities of image
        intensityRange(1) = min(img(:));
        intensityRange(2) = max(img(:));
        % output possible range 
        fprintf('Channel %i. Range in third coordinate: [%i %i].\n',iChan,intensityRange(1),intensityRange(2));
        % request input from user
        userRange = input('Please provide range of intensities to image: ');
        while userRange(1)>userRange(2) 
            userRange = input('The maximum cannot be smaller than the minimum. Please try again: ');
        end
        irange(iChan,:) = [userRange(1) userRange(2)];
        fprintf('\n');
        
    else
        irange(iChan,:) = opts.contrast(iChan,:);
    end
    rgbImg(:,:,chanOrder(iChan)) = imadjust(img, irange(iChan,:), []);

end

% do chromatic shift correction of image
for iChan = opts.imageChans(opts.imageChans~=coordSysChan)
    
    [ rgbImgShift(:,:,chanOrder(coordSysChan)) , rgbImgShift(:,:,chanOrder(iChan)) , ~ ] = ...
        chrsComputeCorrectedImage(rgbImg(:,:,chanOrder(coordSysChan)), rgbImg(:,:,chanOrder(iChan)),...
        chrShift{coordSysChan,iChan},'coords',opts.coords,'pixelSize',pixelSize,...
        'transpose',opts.transpose,'subpixelate',opts.subpixelate);
end

%% Plotting the image

if ~opts.withinFig
    figure(1)
    clf
end
if length(opts.imageChans) == 1
    imshow(rgbImg(:,:,chanOrder(opts.imageChans)));
else
    imshow(rgbImgShift)
end
if opts.crop == -1
    hold on
    rectangle('Position',job.ROI(jobID).crop*opts.subpixelate,'EdgeColor','w')
end

%% Overlays
indent = 10;

% scalebar
barLength = opts.scaleBarSize/pixelSize(1);
labelGap = 6;
barLength = barLength*opts.subpixelate;
indent = indent*opts.subpixelate;
labelGap = labelGap*opts.subpixelate;
imageSize = cropSize*opts.subpixelate;

if barLength > 0
    line([indent indent+barLength], [imageSize(1)-indent imageSize(1)-indent],...
      'Color','w','LineWidth',5);
    if opts.scaleBarLabel
        scaleLabel = [num2str(opts.scaleBarSize) ' µm'];
        text(indent+(barLength/2),imageSize(1)-indent-labelGap,scaleLabel,...
          'Color','w','HorizontalAlignment','center','FontSize',20);
    end
end

% text in each other corner
if ~isempty(opts.textNorthWest)
    % print label requested in top right
    if iscell(opts.textNorthWest)
        vertIndent = indent+((length(opts.textNorthWest)-1)*5);
    else
        vertIndent = indent;
    end
    text(indent,vertIndent,opts.textNorthWest,...
        'Color','w','FontSize',25,'HorizontalAlignment','left');
end
if ~isempty(opts.textNorthEast)
    % print label requested in top right
    if iscell(opts.textNorthEast)
        vertIndent = indent+((length(opts.textNorthEast)-1)*5);
    else
        vertIndent = indent;
    end
    text(imageSize(2)-indent,vertIndent,opts.textNorthEast,...
        'Color','w','FontSize',25,'HorizontalAlignment','right');
end
if ~isempty(opts.textSouthEast)
    % print label requested in top right
    if iscell(opts.textSouthEast)
        vertIndent = imageSize(1)-indent-((length(opts.textSouthEast)-1)*5);
    else
        vertIndent = imageSize(1)-indent;
    end
    text(imageSize(2)-indent,vertIndent,opts.textSouthEast,...
        'Color','w','FontSize',25,'HorizontalAlignment','right');
end

end