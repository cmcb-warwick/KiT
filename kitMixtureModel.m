function job = kitMixtureModel(job,movie,localMaxima,channel)
% Refine spots in 3D using Gaussian mixture model fits
%
% Created by: Jacques Boisvert and K. Jaqaman
% Modified by: J. W. Armond and C. A. Smith
% Copyright (c) 2016 C. A. Smith

%% Input + initialization

dataStruct = job.dataStruct{channel};
options = job.options;
dataStruct.failed = 0;

% get pixel and chromatic shift information
pixelSize = job.metadata.pixelSize;
chrShift = job.options.chrShift.result{options.coordSystemChannel,channel};
% calculate alphaA based on pixel size
% (standards: 0.05 at 138.1nm; 0.5 at 69.4nm)
options.mmf.alphaA = options.mmf.alphaA.*(options.mmf.alphaA + (0.1381-pixelSize(1))*8)/0.05;
% NOTE TO SELF: Use scaling factor = 6.5 for 0.5, or 8 for 0.6

% get number of frames
nFrames = job.metadata.nFrames;
is3D = job.metadata.is3D;
ndims = 2 + is3D;
%get initial guess of PSF sigma
filterPrm = dataStruct.dataProperties.FILTERPRM;
%initialize some variables
emptyFrames = [];

% get neighbour information
neighbour = strcmp(options.spotMode{channel},'neighbour');
if neighbour
  if isempty(options.neighbourSpots.timePoints{channel})
    frames = 1:nFrames;
  else
    frames = options.neighbourSpots.timePoints{channel};
  end
else
  frames = 1:nFrames;
end
% Go over all frames and register empty frames.
for iImage = frames
  %if there are no cands, register that this is an empty frame
  if isempty(localMaxima(iImage).cands)
    emptyFrames = [emptyFrames iImage]; %#ok<AGROW>
  end
end

%make a list of images that have local maxima, ensuring it is a row
goodImages = setxor(frames,emptyFrames);
goodImages = goodImages(:)';

% get psf sigma from filterPrm
if is3D
  psfSigma = [filterPrm(1) filterPrm(3)];
else
  psfSigma = filterPrm(1);
end

%initialize initCoord
initCoord(frames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0,'amp',[],'bg',[]);

kitLog('Refining particles using mixture-model fitting');
prog = kitProgress(0);

%if neighbour, firstly go over all empty images to give nan structures
if neighbour
  for iImage = emptyFrames
    % Get number of spots in coordinate system channel
    nSpotsOrig = job.dataStruct{options.coordSystemChannel}.initCoord(iImage).nSpots;
    % Create an array of nans and insert neighbour coordinates per ID
    initCoord(iImage).allCoord = nan(nSpotsOrig,2*ndims);
    initCoord(iImage).allCoordPix = nan(nSpotsOrig,2*ndims);
    initCoord(iImage).amp = nan(nSpotsOrig,3);
    initCoord(iImage).bg = nan(nSpotsOrig,2);
  end
end

% Go over all non-empty images ...
for iImage = goodImages
  % Get frame.
  imageRaw = movie(:,:,:,iImage);
  if neighbour
    % Get number of spots in coordinate system channel
    nSpotsOrig = job.dataStruct{options.coordSystemChannel}.initCoord(iImage).nSpots;
  end
  % Get candidate maxima, and append spotIDs to the end
  cands = localMaxima(iImage).cands;
  if neighbour
    cands(:,end+1) = localMaxima(iImage).spotID;
  else
    cands(:,end+1) = NaN;
  end

  % Fit with mixture-models.
  startTime = clock;
  
  % Impose opts with channel-specific mmfAddSpots and alphas.
  opts = options;
  opts.mmf.addSpots = opts.mmf.addSpots*~strcmp(options.spotMode{channel},'neighbour');
  opts.mmf.alphaA = opts.mmf.alphaA(channel);
  opts.mmf.alphaD = opts.mmf.alphaD(channel);
  opts.mmf.alphaF = opts.mmf.alphaF(channel);
  % Run mixture model fitting.
  [coordList,spotID,ampList,bgList,rejects] = mixtureModelFit(cands,imageRaw,psfSigma,opts);
  
  elapsedTime = etime(clock,startTime);
  if options.mmf.maxMmfTime > 0 && elapsedTime > options.mmf.maxMmfTime
    warning('Mixture-model fitting taking excessive time (%g min per frame). Aborting.',elapsedTime/60);
    dataStruct.failed = 1;
    break;
  end
  nSpots = size(coordList,1);
  if ~is3D
    coordList = [coordList(:,1:2) zeros(nSpots,1) coordList(:,3:4) zeros(nSpots,1)];
  end

  % Visualize final result.
  if options.debug.showMmfFinal ~= 0
    % If 3D image, max project.
    img = max(imageRaw,[],3);
    figure(1);
    imshow(img,[]);

    % Plot image and overlay spots.
    hold on;
    plot(cands(:,2),cands(:,1),'b+');
    if ~isempty(rejects.amp)
      plot(rejects.amp(:,1),rejects.amp(:,2),'gx');
    end
    if ~isempty(rejects.dist)
      plot(rejects.dist(:,1),rejects.dist(:,2),'yx');
    end
    if ~isempty(coordList)
      plot(coordList(:,1),coordList(:,2),'rx');
    end

    title('MMF fitted particles (r), cands (b), amp rej (g), dist rej (y)');
    hold off;
    drawnow;
    switch options.debug.showMmfFinal
      case -1
        pause;
      case -2
        keyboard;
    end
  end

  if options.debug.showMmfPvals ~= 0
    figure(2);
    subplot(2,1,1);
    if ~isempty(rejects.amp)
      histogram(rejects.amp(:,ndims+1));
    end
    title('Amplitude reject p-vals');
    subplot(2,1,2);
    if ~isempty(rejects.dist)
      histogram(rejects.dist(:,ndims));
    end
    title('Distance reject p-vals');
    drawnow;
    switch options.debug.showMmfPvals
      case -1
        pause;
      case -2
        keyboard;
    end
  end
  
  %save results
  initCoord(iImage).nSpots = nSpots;
  if neighbour
    % Create an array of nans and insert neighbour coordinates per ID
    initCoord(iImage).allCoord = nan(nSpotsOrig,2*ndims);
    initCoord(iImage).allCoordPix = nan(nSpotsOrig,2*ndims);
    initCoord(iImage).amp = nan(nSpotsOrig,3);
    initCoord(iImage).bg = nan(nSpotsOrig,2);
    for iSpot = 1:nSpots
      initCoord(iImage).allCoord(spotID(iSpot),:) = (coordList(iSpot,:).*repmat(pixelSize,1,2)) - [chrShift(1:ndims) 0 0 0];
      initCoord(iImage).allCoordPix(spotID(iSpot),:) = coordList(iSpot,:);
      initCoord(iImage).amp(spotID(iSpot),:) = ampList(iSpot,:);
      initCoord(iImage).bg(spotID(iSpot),:) = bgList(iSpot,:);
    end
  
  else
      
    initCoord(iImage).allCoordPix = coordList;
    initCoord(iImage).amp = ampList;
    initCoord(iImage).bg = bgList;

    %check whether frame is empty
    if initCoord(iImage).nSpots == 0
      emptyFrames = [emptyFrames iImage]; %#ok<AGROW>
      initCoord(iImage).allCoord = initCoord(iImage).allCoordPix;
    else
      % calc real space coordinates and correct for chromatic shift
      initCoord(iImage).allCoord = initCoord(iImage).allCoordPix .* repmat(pixelSize,initCoord(iImage).nSpots,2);
      initCoord(iImage).allCoord(:,1:ndims) = initCoord(iImage).allCoord(:,1:ndims) - repmat(chrShift(1:ndims),initCoord(iImage).nSpots,1);
    end
  
  end

  %display progress
  prog = kitProgress(iImage/length(goodImages),prog);
end

%% Post-processing

%sort list of empty frames, keep only unique frames
emptyFrames = unique(emptyFrames);

%store empty frames and frames where detection failed in structure
%exceptions
exceptions = struct('emptyFrames',emptyFrames);

% save results
initCoord(1).exceptions = exceptions;
initCoord(1).localMaxima = localMaxima;
dataStruct.dataProperties.psfSigma = psfSigma;
dataStruct.initCoord = initCoord;

dataStruct.failed = dataStruct.failed || length(emptyFrames) == nFrames;

job.dataStruct{channel} = dataStruct;
