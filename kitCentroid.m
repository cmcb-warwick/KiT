function job=kitCentroid(job,movie,localMaxima,channel)
% Refine spots in 3D using centroids
%
% Copyright (c) 2015 Jonathan W. Armond

dataStruct = job.dataStruct{channel};
ndims = 2 + (size(movie,3)==1);
filters = createFilters(ndims,dataStruct.dataProperties);
nFrames = length(localMaxima);
opts = job.options;

% get pixel and chromatic shift information
pixelSize = job.metadata.pixelSize;
chrShift = job.options.chrShift.result{job.options.coordSystemChannel,channel};

% find whether this is a neighbour-localised channel or not
neighbour = strcmp(job.options.spotMode{channel},'neighbour');

% get halfPatchSize to adjust centroid result. The center pixel of a 5x5x5
% array is (3,3,3), thus, we have to subtract that from centroid coordinates
halfPatchSize = dataStruct.dataProperties.FILTERPRM(4:6)/2+0.5;
kitLog('Refining spots using centroids');
prog = kitProgress(0);
for t = 1:nFrames
  
  if neighbour
    % get spotIDs from localMaxima structure, and number of
    % coordSystemChannel spots
    spotIDs{t} = localMaxima(t).spotID;
    nSpotsOrig(t) = job.dataStruct{job.options.coordSystemChannel}.initCoord(t).nSpots;
  end  
    
  for iSpot = 1:size(localMaxima(t).cands,1)
    % read volume around coordinate
    patch = stamp3d(movie(:,:,:,t),filters.signalP,localMaxima(t).cands(iSpot,:),0);

    %HLE,KJ - calculate low-index edge patch adjustment if relevant
    edgeAdjustTmp = localMaxima(t).cands(iSpot,:) - halfPatchSize;
    edgeAdjustTmp = abs(edgeAdjustTmp) .* (edgeAdjustTmp<0);

    % subpixel coord is integer coord plus centroid (subtract
    % halfPatchSize so that center coordinate of patch is (0,0,0))
    centroid = centroid3D(patch);
    if any(centroid>2*halfPatchSize | centroid < 0)
      % Try again with exponent 2.
      centroid = centroid3D(patch,2);
      if any(centroid>2*halfPatchSize | centroid < 0)
        warning('Centroid fitting error: initCoordTmp(%d) t=%d',iSpot,t);
      end
    end
    localMaxima(t).cands(iSpot,:) = localMaxima(t).cands(iSpot,:) + ...
        centroid - halfPatchSize;

    %HLE,KJ - correct position due to patch truncation due to proximity
    %to a low-index edge
    localMaxima(t).cands(iSpot,:) = localMaxima(t).cands(iSpot,:) + ...
        edgeAdjustTmp;
    
    % chromatic shift correction
    localMaxima(t).cands(iSpot,:) = localMaxima(t).cands(iSpot,:) - ...
        chrShift(1:3);

    % amplitude guess is integral.
    localMaxima(t).candsAmp(iSpot) = nanmean(patch(:));
  end

  % display progress
  prog = kitProgress(t/nFrames,prog);
end

% Save results
initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0,'amp',[],'bg',[]);
initCoord(1).localMaxima = localMaxima;

% loop and store only good locMax.
for t=1:nFrames
  
  % Count spots
  initCoord(t).nSpots = size(localMaxima(t).cands,1);

  if neighbour
    
    initCoord(t).allCoordPix = nan(nSpotsOrig(t),6);
    % Store pixel coords in image coordinate system (i.e. x == cols, y ==
    % rows). Uncertainty is 0.25 pix.
    initCoord(t).allCoordPix(spotIDs{t},:) = [localMaxima(t).cands(:,[2 1 3]) ...
       0.25*ones(initCoord(t).nSpots,3)];
  
    % store estimated amplitude and noise
    %initCoord(t).initAmp = initCoordRaw{t}(spotIDs(t),4:5);
    
    initCoord(t).amp = nan(nSpotsOrig(t),2);
    % store integral amplitude
    initCoord(t).amp(spotIDs{t},:) = [localMaxima(t).candsAmp zeros(initCoord(t).nSpots,1)];
    
    % store coords in microns and correct
    initCoord(t).allCoord = bsxfun(@times,initCoord(t).allCoordPix,...
      repmat(pixelSize,[1 2]));
  
  else
      
    % Store pixel coords in image coordinate system (i.e. x == cols, y ==
    % rows). Uncertainty is 0.25 pix.  
    initCoord(t).allCoordPix = [localMaxima(t).cands(:,[2 1 3]) ...
       0.25*ones(initCoord(t).nSpots,3)];
   
    % store estimated amplitude and noise
    %initCoord(t).initAmp = initCoordRaw{t}(goodIdxL,4:5);
    
    % store integral amplitude
    initCoord(t).amp = [localMaxima(t).candsAmp zeros(initCoord(t).nSpots,1)];

    % store coords in microns and correct
    initCoord(t).allCoord = bsxfun(@times,initCoord(t).allCoordPix,...
      repmat(pixelSize,[1 2]));
    
  end

  % Visualize final result.
  if opts.debug.showCentroidFinal ~= 0
    % If 3D image, max project.
    raw = kitReadImageStack(reader,metadata,t,channel,job.ROI.crop);
    img = max(raw,[],3);
    figure(1);
    imshow(img,[]);

    % Plot image and overlay spots.
    hold on;
    if ~isempty(localMaxima(t).cands)
      plot(localMaxima(t).cands(:,2),localMaxima(t).cands(:,1),'b+');
    end
    if initCoords(t).nSpots > 0
      plot(initCoord(t).allCoordPix(:,1), ...
           initCoord(t).allCoordPix(:,2),'rx');
    end

    title('Centroid fitted spots (r), cands (b)');
    hold off;
    drawnow;
    switch opts.debug.showCentroidFinal
      case -1
        pause;
      case -2
        keyboard;
    end
  end
end

dataStruct.initCoord = initCoord;
dataStruct.failed = 0;
job.dataStruct{channel} = dataStruct;

