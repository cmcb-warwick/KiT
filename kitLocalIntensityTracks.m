function job=kitLocalIntensityTracks(job, reader, metadata, channel, varargin)
% KITLOCALINTENSITY Measure intensity local to kinetochores

% Set defaults.
opts.maskRadius = job.options.maskRadius;
opts.maskShape = job.options.maskShape;
opts.maskConeAngle = job.options.maskConeAngle;
opts.otherSpotSearchRadius = job.options.otherSpotSearchRadius;
opts.photobleachCorrect = job.options.photobleachCorrect;
opts.poleShift = job.options.poleShift;
opts.gaussFilterSpots = job.options.gaussFilterSpots;

% Process options.
opts = processOptions(opts, varargin{:});

initCoord = job.dataStruct{channel}.initCoord;
trackList = job.dataStruct{channel}.trackList;
planeFit = job.dataStruct{channel}.planeFit;
nTracks = length(trackList);
poleShiftPixels = round(opts.poleShift / metadata.pixelSize(1));

nFrames = metadata.nFrames;
nChannels = metadata.nChannels;
trackInt(1:nTracks) = struct(...
    'intensity',nan(nFrames,nChannels),...
    'intensity_median',nan(nFrames,nChannels),...
    'intensity_min',nan(nFrames,nChannels),...
    'intensity_max',nan(nFrames,nChannels),...
    'intensity_ratio',nan(nFrames,nChannels),...
    'maskCoord',nan(nFrames,3),...
    'maxCoord',nan(nFrames,3*nChannels),...
    'angleToMax',nan(nFrames,nChannels),...
    'distToMax',nan(nFrames,nChannels));


% Create mask.
r = round(opts.maskRadius / metadata.pixelSize(1));
if ~strcmp(opts.maskShape, 'cone')
  se = strel('disk', r, 0);
  asymMask = 0;
  mask = double(se.getnhood());
  mask(mask == 0) = nan;
  if strcmp(opts.maskShape, 'semicircle')
    mask(r+2:end,:) = nan;
    asymMask = 1;
  end
else
  mask = ones(r);
  % Make conical mask.
  ycut = tand(opts.maskConeAngle)*(0:r-1);
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

maskWarning=0;

% Gaussian filter.
hgauss = fspecial('gauss');

% Read frame by frame.
prog = kitProgress(0);
for t=1:nFrames
  for c=1:nChannels
    stack = kitReadImageStack(reader, metadata, t, c, job.crop, 0);
    if opts.gaussFilterSpots
      % Gaussian filter each z-plane.
      for z=1:size(stack,3)
        stack(:,:,z) = imfilter(stack(:,:,z), hgauss, 'replicate');
      end
    end

    if t==1
      % Take overall fluorescence estimate.
      back(c) = mean(stack(:));

      % Estimate background and signal.
      [backMode(c),sigVal,imgHist] = estBackgroundAndSignal(stack);
      intensityDistF(:,c) = imgHist(:,2);
      intensityDistX(:,c) = imgHist(:,1);

      % Estimate background as difference between signal mode and background mode.
      backDiff(c) = sigVal - backMode(c);
    end

    for j=1:nTracks
      % Map track back to coords. FIXME use trackList instead
      pixIdx = trackList(j).featIndx(t);

      if pixIdx > 0 && pixIdx < size(initCoord(t).allCoordPix,1)
        % Read off intensity in radius r around spot.
        pixCoords = initCoord(t).allCoordPix(pixIdx,1:3);
        if opts.otherSpotSearchRadius > 0 && c ~= channel
          % Try to locate nearby spot in other channel to use as centre
          % of radius to read intensity.
          otherCoordPix = job.dataStruct{c}.initCoord(t).allCoordPix(:,1:3);
          otherCoord = job.dataStruct{c}.initCoord(t).allCoord(:,1:3);
          coords = initCoord(t).allCoord(pixIdx,1:3);
          dists = sqrt(sum((otherCoord-repmat(coords,size(otherCoord,1),1)).^2,2));
          [nearDist, nearIdx] = min(dists);
          if nearDist > opts.otherSpotSearchRadius
            continue
          end
          % Use the near other coords to read pixel
          % intensity.
          pixCoords = otherCoordPix(nearIdx, :);
        end

        if sum(isnan(pixCoords)) == 0
          if trackList(j).attach ~= 0 && ...
                size(pixCoords,2)==size(planeFit(t).planeVectors,2)
            % Shift mask toward pole, if known sister.
            maskCoords = pixCoords + [sign(trackList(j).attach),0,0] * ...
                poleShiftPixels*planeFit(t).planeVectors;
          else
            maskCoords = pixCoords;
          end

          % Extract intensity.
          x = max(1,round(maskCoords(1)));
          y = max(1,round(maskCoords(2)));
          z = max(1,round(maskCoords(3)));
          if z > size(stack,3)
            warning('Spot outside frame boundaries');
            continue;
          end
          if size(stack,3)>1
            imgPlane = stack(:,:,z);
          else
            imgPlane = stack;
          end
          [my,mx] = size(imgPlane);
          % If too close to edge, skip.
          if x <= r || y <= r || x >= mx-r || y >= my-r
            continue;
          end

          % If asymmetric mask, choose left/right mask based on pole
          % attachment.
          if asymMask
            if trackList(j).attach > 0
              mask = rMask;
            else
              mask = lMask;
            end
            if ~isempty(planeFit(t).planeVectors)
              % Rotate mask into coordinate system.
              angle = acos(planeFit(t).planeVectors(1,1))*180/pi;
              mask = imrotate(mask, angle, 'nearest', 'crop');
            else
              if ~maskWarning
                warning('No rotation vector for mask');
                maskWarning=1;
              end
            end
          end
          % Coordinates are in image system, so plot(x,y) should draw the
          % spots in the correct place over the image. However, the image
          % matrix is indexed by (row,col) => (y,x).
          maskImg = mask .* imgPlane(y-r:y+r, x-r:x+r);
          nonNanPix = maskImg(~isnan(maskImg));
          trackInt(j).intensity(t,c) = mean(nonNanPix);
          trackInt(j).intensity_median(t,c) = median(nonNanPix);
          trackInt(j).intensity_max(t,c) = max(nonNanPix);
          trackInt(j).intensity_min(t,c) = min(nonNanPix);
          trackInt(j).intensity_ratio(t,c) = trackInt(j).intensity_max(t,c) ...
              / trackInt(j).intensity_min(t,c);
          trackInt(j).maskCoord(t,:) = maskCoords;
          %[maxX,maxY] = ind2sub(size(maskImg),maxIdx);
          [~,maxIdx] = max(maskImg(:));
          [maxY,maxX] = ind2sub(size(maskImg),maxIdx);
          trackInt(j).maxCoord(t,3*(c-1)+1:3*c) = [maxX+x-r-1,maxY+y-r-1,z];

          % Calculate angle between spot and max point.
          vector = pixCoords-trackInt(j).maxCoord(t,3*(c-1)+1:3*c);
          trackInt(j).distToMax(t,c) = norm(vector(1:2),2)*metadata.pixelSize(1);
          % Rotate angle into coordinate system.
          if size(pixCoords,2)==size(planeFit(t).planeVectors,2)
            vector = vector*planeFit(t).planeVectors;
          end
          angle = -atan2(vector(2),-vector(1)); % See doc atan2 for diagram.
          trackInt(j).angleToMax(t,c) = angle;
        end
      end
    end
  end

  % Report progress.
  prog = kitProgress(t/nFrames, prog);
end

if ~isempty(opts.photobleachCorrect)
  pbProfile = opts.photobleachCorrect;
  % Correct intensities using photobleching profile.
  if isscalar(pbProfile)
    % Compute photobleach from entire image.
    pbProfile=zeros(job.metadata.nFrames,job.metadata.nChannels);
    for c=1:metadata.nChannels
      pbProfile(:,c) = kitIntensityDistn(job,reader,metadata,c,[],[],1,job.crop);
    end
  elseif size(pbProfile,2) == 1
    % If not invidual channels specified, use the same.
    pbProfile = repmat(pbProfile,metadata.nChannels,1);
  end
  t=((1:size(pbProfile,1))-1)';

  for c=1:metadata.nChannels
    % Normalize PB if necessary.
    if pbProfile(1,c) ~= 1 || ~all(pbProfile(:,c) < 1)
      pbProfile(:,c) = pbProfile(:,c)/pbProfile(1,c);
    end

    if numel(t) > 4 && license('test','Curve_Fitting_Toolbox')
      % Fit double exp.
      pbF = fit(t,pbProfile(:,c),'exp2');
    else
      kitLog('Warning: Not correcting for photobleach');
      break;
    end

    % Correct.
    for j=1:nTracks
      trackInt(j).intensity(:,c) = trackInt(j).intensity(:,c)./pbF(t);
      trackInt(j).intensity_median(:,c) = trackInt(j).intensity_median(:,c)./pbF(t);
      trackInt(j).intensity_max(:,c) = trackInt(j).intensity_max(:,c)./pbF(t);
      trackInt(j).intensity_min(:,c) = trackInt(j).intensity_min(:,c)./pbF(t);
    end
  end
end

% record background intensity
cellInt.back = back;
cellInt.backMode = backMode;
cellInt.backDiff = backDiff;
cellInt.maskPixelRadius = r;
cellInt.intensityDistF = intensityDistF;
cellInt.intensityDistX = intensityDistX;

% return data
job.dataStruct{channel}.cellInt = cellInt;
job.dataStruct{channel}.trackInt = trackInt;

reader.close();
