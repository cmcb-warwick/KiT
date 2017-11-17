function job = kitMixtureModel(job,movie,localMaxima,channel)
% Refine spots in 3D using Gaussian mixture model fits
%
% Created by: Jacques Boisvert and K. Jaqaman
% Modified by: J. W. Armond
% Copyright (c) 2014 Jonathan W. Armond

%% Input + initialization

dataStruct = job.dataStruct{channel};
options = job.options;
dataStruct.failed = 0;

pixelSize = job.metadata.pixelSize;
% calculate alphaA based on pixel size
% (standards: 0.05 at 138.1nm; 0.5 at 69.4nm)
options.alphaA = options.alphaA.*(options.alphaA + (0.1381-pixelSize(1))*8)/0.05;
% NOTE TO SELF: Use scaling factor = 6.5 for 0.5, or 8 for 0.6

% get number of frames
nFrames = job.metadata.nFrames;
is3D = job.metadata.is3D;
ndims = 2 + is3D;
%get initial guess of PSF sigma
filterPrm = dataStruct.dataProperties.FILTERPRM;
%initialize some variables
emptyFrames = [];

[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);

% Go over all frames and register empty frames.
for iImage = 1 : nFrames
  %if there are no cands, register that this is an empty frame
  if isempty(localMaxima(iImage).cands)
    emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
  end
end

%make a list of images that have local maxima
goodImages = setxor(1:nFrames,emptyFrames,'legacy');

% get psf sigma from filterPrm
if is3D
  psfSigma = [filterPrm(1) filterPrm(3)];
else
  psfSigma = filterPrm(1);
end

%initialize initCoord
initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0,'amp',[],'bg',[]);

kitLog('Refining particles using mixture-model fitting');
prog = kitProgress(0);
%go over all non-empty images ...
for iImage = goodImages
  % Get frame.
  imageRaw = movie(:,:,:,iImage);

  % Get candidate maxima.
  cands = localMaxima(iImage).cands;
  
  % Impose opts with channel-specific mmfAddSpots and alphas.
  opts = options;
  opts.alphaA = opts.alphaA(channel);
  opts.alphaD = opts.alphaD(channel);
  opts.alphaF = opts.alphaF(channel);

  % Fit with mixture-models.
  startTime = clock;
  [coordList,ampList,bgList,rejects] = mixtureModelFit(cands,imageRaw,psfSigma,opts);
  elapsedTime = etime(clock,startTime);
  if options.maxMmfTime > 0 && elapsedTime > options.maxMmfTime
    warning('Mixture-model fitting taking excessive time (%g min per frame). Aborting.',elapsedTime/60);
    dataStruct.failed = 1;
    break;
  end
  nSpots = size(coordList,1);
  if ~isempty(coordList) && ~is3D
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
  initCoord(iImage).allCoordPix = coordList;
  initCoord(iImage).amp = ampList;
  initCoord(iImage).bg = bgList;


  %check whether frame is empty
  if initCoord(iImage).nSpots == 0
    emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
    initCoord(iImage).allCoord = initCoord(iImage).allCoordPix;
  else
    % calc real space coordinates
    initCoord(iImage).allCoord = initCoord(iImage).allCoordPix .* repmat(job.metadata.pixelSize,initCoord(iImage).nSpots,2);
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
