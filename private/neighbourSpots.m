function [spots,spotIDs] = neighbourSpots(movie,refDataStruct,channel,metadata,options)
% NEIGHBOURSPOTS Finds spots in 3D in neighbouring channels using Gaussian
% mixture model fitting
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2017 C. A. Smith

%% Input + initialization

% get all required coordinate information
refInitCoord = refDataStruct.initCoord;

% get pixel and chromatic shift information
pixelSize = metadata.pixelSize(1:3);
chrShift = options.chrShift.result{options.coordSystemChannel,channel}(:,1:3);
chrShift = chrShift./pixelSize; chrShift(1:2) = chrShift([2 1]);
% image size
imageSize = size(movie);
nFrames = metadata.nFrames;
nPlanes = imageSize(3);
% get channel orientation (1 if channel outside reference, 0 otherwise)
refChanPos = find(options.neighbourSpots.channelOrientation==options.coordSystemChannel);
chanPos = find(options.neighbourSpots.channelOrientation==channel);
chanOrient = (chanPos>refChanPos);

% get list of frames over which to look for neighbours
if isempty(options.neighbourSpots.timePoints{channel})
    timePoints = 1:nFrames;
else
    timePoints = options.neighbourSpots.timePoints{channel};
end
if isfield(refDataStruct.initCoord(1),'exceptions')
    emptyFrames = refDataStruct.initCoord(1).exceptions.emptyFrames;
else
    emptyFrames = [];
end
emptyFrames = union(emptyFrames, setxor(timePoints,1:nFrames));
timePoints = setxor(emptyFrames,1:nFrames); % final list of frames
timePoints = timePoints(:)'; % ensure is a row vector - setxor has undergone a change over MATLAB versions

%get number of z-slices
if isempty(options.neighbourSpots.zSlices{channel})
    zSlices = 1:nPlanes;
else
    zSlices = options.neighbourSpots.zSlices{channel};
end

%turn warnings off
warningState = warning();

%% Get frame with plane fits

if ~strcmp(options.jobProcess,'chrshift')
  % get plane fit information from reference channel
  planeFit = refDataStruct.planeFit;
 
  % find frames with and without plane fit
  framesNoPlane = [];
  for iFrame = 1:nFrames
    if isempty(planeFit(iFrame).plane)
      framesNoPlane = [framesNoPlane iFrame];
    end
  end
else
  framesNoPlane = 1;
end

%% Create Mask

r = round(options.neighbourSpots.maskRadius/pixelSize(1));
if ~strcmp(options.neighbourSpots.maskShape, 'cone')
  se = strel('disk', r, 0);
  asymMask = 0;
  mask = double(se.getnhood());
  mask(mask == 0) = nan;
  if strcmp(options.neighbourSpots.maskShape, 'semicirc')
    mask(r+2:end,:) = nan;
    asymMask = 1;
  end
else
  mask = ones(r);
  % Make conical mask.
  ycut = tand(options.neighbourSpots.maskConeAngle)*(0:r-1);
  for i=1:r
    ycutFlr = floor(ycut(i));
    if ycutFlr>=1
      mask(1:ycutFlr,i) = nan;
    end
  end
  % Add asymmetry and orient.
  mask = [flipud(mask); nan(r-1,r)];
  % Mirror
  mask = [fliplr(mask) mask(:,2:r)];
  asymMask = 1;
end

if asymMask == 1
  lMask = mask; % For -ve x pole attached.
  rMask = flipud(mask); % For +ve x pole attached.
end
if ~isempty(framesNoPlane)
    se = strel('disk', r, 0);
    noPlaneMask = double(se.getnhood());
    noPlaneMask(noPlaneMask == 0) = nan;
end

maskWarning=0;

%% Local maxima detection relative to reference channel, based on job process

switch options.jobProcess
    
  case 'zandt'
      
    % get trackList and sisterList from reference channel
    refTrackList = refDataStruct.trackList;
    refSisterList = refDataStruct.sisterList;

    % get number of sisters
    if ~isempty(refSisterList(1).trackPairs)
        nSisters = length(refSisterList);
    else
        warning('No sisterList found for channel %d. Tracking failed.',options.coordSystemChannel)
        spots = repmat({[]},nFrames,1);
        spotIDs = repmat({[]},nFrames,1);
        return
    end
    
    % create structure to store sister number and track numbers numbers
    referenceIDs = nan(nSisters,3);
    % create spots and spotIDs structures
    spots = repmat({[]},nFrames,1);
    spotIDs = repmat({[]},nFrames,1);

    for iSisPair = 1:nSisters

        % get sisterList, and its track and spot IDs
        sisList = refSisterList(iSisPair);
        iTracks = refSisterList(1).trackPairs(iSisPair,1:2);
        spotIDbySister{iSisPair}  = [refTrackList(iTracks(1)).featIndx refTrackList(iTracks(2)).featIndx];
        referenceIDs(iSisPair,:) = [iSisPair iTracks];

        for iFrame = timePoints
          % get initCoord for this timepoint
          initCoord = refInitCoord(iFrame);

          for iSis = 1:2

              trackID = iTracks(iSis);
              spotID = spotIDbySister{iSisPair}(iFrame,iSis);
              % check that this coordinate has a spotID
              if isnan(spotID)  
                continue
              end

              % get coordinate information for this each sister, then transpose
              % (the results of MMF are transposed, so need to be transposed back)
              coords = initCoord.allCoordPix(spotID,[2 1 3]);
              
              % chromatic shift coordinates to new channel's frame, then
              % find nearest whole pixel
              coords = coords + chrShift;
              coords = round(coords);
              
              % define and check the image range
              range = [coords(1)-r coords(1)+r;...
                       coords(2)-r coords(2)-r;...
                       coords(3)   coords(3) ];

              if any(isnan(range(:))) || any(range(:,1)<=0) || any(range(:,2)>imageSize')
                  continue
              elseif range(3,1)<min(zSlices) || range(3,2)>max(zSlices)
                  continue
              end
                  
              % realign mask into correct orientation for this track
              if sum(framesNoPlane == iFrame) > 0
                  mask = noPlaneMask;
              else
                  if refTrackList(trackID).attach > 0
                    if chanOrient
                      mask = rMask;
                    else
                      mask = lMask;
                    end
                  elseif ~chanOrient
                      mask = rMask;
                  else
                      mask = lMask;
                  end
                  
                  if ~isempty(planeFit(iFrame).planeVectors)
                    % rotate mask into coordinate system
                    angle = acos(planeFit(iFrame).planeVectors(1,1))*180/pi;
                    mask = imrotate(mask, angle, 'nearest', 'crop');
                  else
                    if ~maskWarning
                      warning('No rotation vector for mask.');
                      maskWarning=1;
                    end
                  end
              end
                  
              % get frame
              image = movie(:,:,:,iFrame);
              
              % get intensity values of mask-derived spot vicinity
              imageMask = mask .* image(range(1,1):range(1,2), ...
                  range(2,1):range(2,2), range(3,1):range(3,2));

              % find local maximum intensity
              locMax1DIndx = find(imageMask==nanmax(imageMask(:)));
              if isempty(locMax1DIndx)
                  continue
              elseif length(locMax1DIndx)>1
                  locMax1DIndx = locMax1DIndx(1); % IDEALLY SHOULD BE THE SPOT CLOSEST TO THE ORIGINAL SPOT
              end

              [locMaxCrd(1),locMaxCrd(2),locMaxCrd(3)] = ind2sub([2*r+1 2*r+1 1],locMax1DIndx);
              % correct coordinates to full image
              locMaxCrd(1) = locMaxCrd(1)+coords(1)-(r+1);
              locMaxCrd(2) = locMaxCrd(2)+coords(2)-(r+1);
              locMaxCrd(3) = coords(3);

              if any(locMaxCrd>imageSize)
                  continue
              end

              % compile coordinates
              spots{iFrame} = [spots{iFrame}; locMaxCrd];
              spotIDs{iFrame} = [spotIDs{iFrame}; spotID];

          end %iSis
        end %iFrame

    end %iSisPair
    
  case {'zonly','chrshift'}
    
    % get number of spots from reference initCoord
    nSpots = refInitCoord.nSpots;
    % create spots and spotIDs structures
    spots = {[]};
    spotIDs = {[]};

    for iSpot = 1:nSpots

      % get coordinate information
      coords = refInitCoord(1).allCoordPix(iSpot,[2 1 3]);
      
      % chromatic shift coordinates to new channel's frame, then
      % find nearest whole pixel
      coords = coords + chrShift;
      coords = round(coords);

      % define and check the image range
      range = [coords(1)-r coords(1)+r;...
               coords(2)-r coords(2)-r;...
               coords(3)   coords(3) ];
           
      if any(isnan(range(:))) || any(range(:,1)<=0) || any(range(:,2)>imageSize')
          continue
      end
  
      % get frame
      image = movie(:,:,:,1);
        
      % get intensity values of mask-derived spot vicinity
      imageMask = mask .* image(range(1,1):range(1,2), ...
          range(2,1):range(2,2), range(3,1):range(3,2));

      % find local maximum intensity
      locMax1DIndx = find(imageMask==nanmax(imageMask(:)));
      if isempty(locMax1DIndx)
          continue
      elseif length(locMax1DIndx)>1
          locMax1DIndx = locMax1DIndx(1); % IDEALLY SHOULD BE THE SPOT CLOSEST TO THE ORIGINAL SPOT
      end
        
      [locMaxCrd(1),locMaxCrd(2),locMaxCrd(3)] = ind2sub([2*r+1 2*r+1 1],locMax1DIndx);
      % correct coordinates to full image
      locMaxCrd(1) = locMaxCrd(1)+coords(1)-(r+1);
      locMaxCrd(2) = locMaxCrd(2)+coords(2)-(r+1);
      locMaxCrd(3) = coords(3);
      
      if any(locMaxCrd>imageSize)
          continue
      end

      % compile coordinates for MMF
      spots{1} = [spots{1}; locMaxCrd];
      spotIDs{1} = [spotIDs{1}; iSpot];

    end %iSpot
    
end

%go back to original warnings state
warning(warningState);
