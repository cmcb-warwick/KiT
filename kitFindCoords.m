function job=kitFindCoords(job, reader, channel)
%KITFINDCOORD Find kinetochore coordinates
%
% SYNOPSIS: job=kitFindCoords(job, raw, channel)
%
% INPUT job: Struct containing tracking job setup options.
%            Requires at least the following fields:
%
%       reader: BioFormats reader.
%
%       channel: Channel to find coords in.
%
% OUTPUT job: as input but with updated values.
%
% Copyright (c) 2012 Jonathan W. Armond

% Method of fixing spot locations: centroid or Gaussian MMF, or none.
method = job.options.coordMode{channel};
% Handle old jobset versions.
if ~isfield(job.options,'spotMode')
  spotMode = 'histcut';
else
  % Method of identify first spot candidates: Histogram cut 'histcut',
  % adaptive threshold 'adaptive', or  multiscale wavelet product 'wavelet'.
  spotMode = job.options.spotMode{channel};
end

% Set up data struct.
options = job.options;

nFrames = job.metadata.nFrames;
is3D = job.metadata.is3D;
ndims = 2 + is3D;
filters = createFilters(ndims,job.dataStruct{channel}.dataProperties);

% Read image
movie = kitReadWholeMovie(reader,job.metadata,channel,job.crop,0,1);
[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);
if options.deconvolve
  kitLog('Deconvolving');
  p=kitProgress(0);
  for i=1:nFrames
    movie(:,:,:,i) = deconvlucy(movie(:,:,:,i),job.psf);
    p=kitProgress(i/nFrames,p);
  end
end

% Initialize output structure
localMaxima = repmat(struct('cands',[]),nFrames,1);

% Find candidate spots.
switch spotMode
  case 'histcut'
    kitLog('Detecting particle candidates using unimodal histogram threshold');
    spots = cell(nFrames,1);
    for i=1:nFrames
      img = movie(:,:,:,i);
      spots{i} = histcutSpots(img,options,job.dataStruct{channel}.dataProperties);
    end

  case 'adaptive'
    kitLog('Detecting particle candidates using adaptive thresholding');
    spots = adaptiveSpots(movie,options.adaptiveLambda,options.debug.showAdaptive);

  case 'wavelet'
    kitLog('Detecting particle candidates using multiscale wavelet product');
    if options.waveletLevelAdapt
     options.waveletLevelThresh = waveletAdapt(movie,options);
     job.options = options;
     kitSaveJob(job); % Record used value.
    end

    for i=1:nFrames
      img = movie(:,:,:,i);
      spots{i} = waveletSpots(img,options);
    end
    
  case 'manual'
    kitLog('Detecting particle candidates using manual detection');
    spots = manualDetection(movie,job.metadata,options);

  otherwise
    error('Unknown particle detector: %s',spotMode);
end

nSpots = zeros(nFrames,1);
for i=1:nFrames
  nSpots(i) = size(spots{i},1);

  % Round spots to nearest pixel and limit to image bounds.
  spots{i} = bsxfun(@min,bsxfun(@max,round(spots{i}),1),[imageSizeX,imageSizeY,imageSizeZ]);

  % Store the cands of the current image
  % TODO this is computed in both spot detectors, just return it.
  img = movie(:,:,:,i);
  if verLessThan('images','9.2')
    background = fastGauss3D(img,filters.backgroundP(1:3),filters.backgroundP(4:6));
  else
    background = imgaussfilt3(img,filters.backgroundP(1:3),'FilterSize',filters.backgroundP(4:6));
  end
  localMaxima(i).cands = spots{i};
  spots1D = sub2ind(size(img),spots{i}(:,1),spots{i}(:,2),spots{i}(:,3));
  localMaxima(i).candsAmp = img(spots1D);
  localMaxima(i).candsBg = background(spots1D);

  % Visualize candidates.
  if options.debug.showMmfCands ~= 0
    showSpots(img,spots{i});
    title(['Local maxima cands n=' num2str(size(spots{i},1))]);
    drawnow;
    switch options.debug.showMmfCands
      case -1
        pause;
      case -2
        keyboard;
    end
  end
end
kitLog('Average particles per frame: %.1f +/- %.1f',mean(nSpots),std(nSpots));

% Refine spot candidates.
switch method
  case 'centroid'
    job = kitCentroid(job,movie,localMaxima,channel);
  case 'gaussian'
    job = kitMixtureModel(job,movie,localMaxima,channel);
  case 'norefine'
    % No refinement. Copy localMaxima to initCoords.
    initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0, ...
                                  'amp',[],'bg',[]);
    initCoord(1).localMaxima = localMaxima;
    for i=1:nFrames
      initCoord(i).nSpots = size(localMaxima(i).cands,1);
      initCoord(i).allCoordPix = [localMaxima(i).cands(:,[2 1 3]) ...
                          0.25*ones(initCoord(i).nSpots,3)];
      initCoord(i).allCoord = bsxfun(@times, initCoord(i).allCoordPix,...
        repmat(job.metadata.pixelSize,[1 2]));
      initCoord(i).amp = [localMaxima(i).candsAmp zeros(initCoord(i).nSpots,1)];
    end
    % Store data.
    job.dataStruct{channel}.initCoord = initCoord;
    job.dataStruct{channel}.failed = 0;
  otherwise
    error(['Unknown coordinate finding mode: ' job.coordMode]);
end

