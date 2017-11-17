function K=kitMakeKymograph(job,pairIdx,varargin)
% KITMAKESPOTMOVIE Renders a kymograph of a sister pair
%
%    KITMAKESPOTMOVIE(JOB,PAIRIDX,...) Renders a kymograph of a sister pair
%    associated with JOB. Supply options as string/value pairs following JOB.
%
%    Options, defaults in {}:-
%
%    outfile: {none} or filename. Filename to write movie to.
%
%    points: {100} or integer. Points to interpolate intensity across (x-axis
%    of kymograph).
%
%    trackChannel: {first tracked channel} or channel number. Tracked channel
%    to get spots from.
%
%    saturate: {[1 100]} or 2 element vector or scalar, between 0 and 100. As a
%    scalar, percentage of pixel values, at either end of histogram, to
%    saturate. Increase to enhance constrast more. Alternatively, supply a
%    vector [LOW HIGH] to specify percentage to saturate at low and high pixels
%    values. NB 2 is equivalent to [1 99].
%
%    channelMap: {[2 1 3]} or a 3-element vector containing integers between
%    1 and 3. Controls which channels are shown in which colors. 1 is red, 2
%    is green and 3 is blue.
%
%    scale: {[1 1]} or 2-element vector. Factor to resize final kymograph by
%    in X and Y.
%
%    scaleMethod: {'nearest'} or 'bilinear' or 'bicubic'. Interpolation
%    method for scaling image up.
%
%    fixture: {'sisters'} or 'plate'. Reference point to fix X against,
%    sister pair centre or metaphase plate.
%
%    margin: {1.0} or real. Multiple of extent of sister travel to extend
%    kymograph width by.
%
%    jet: {0} or channel number. Enable jet colorscheme for channel number.
%
%    rotate: {0} or 1. Draw kymograph horizontally.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<2
  error('Must supply JOB and PAIRIDX');
end

% Set defaults.
opts.trackChannel = 0; % Default to first tracked channel.
opts.outfile = [];
opts.saturate = [1 100];
opts.channelMap = [2 1 3]; % Green, red, blue
opts.points = 100;
opts.scale = [1 1];
opts.scaleMethod = 'nearest';
opts.fixture = 'sisters';
opts.margin = 1.0;
opts.jet = 0;
opts.trackIdx = 0; % Use pairIdx as a trackList indices not sisterList.
opts.scaleBar = [1,10]; % [um,sec] for the scale bar
opts.rotate = 0;

% Process options.
opts = processOptions(opts, varargin{:});

if opts.trackChannel == 0
  % Plot first non-empty tracking channel
  opts.trackChannel = find(~cellfun(@isempty,job.dataStruct),1);
end

mapChans = opts.channelMap;

nFrames = job.metadata.nFrames;
nChannels = job.metadata.nChannels;

if size(opts.saturate,1)<nChannels
  opts.saturate=repmat(opts.saturate,[nChannels 1]);
end

% Open movie.
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie));

ds = job.dataStruct{opts.trackChannel};
if ~opts.trackIdx && (pairIdx > length(ds.sisterList) || pairIdx < 1)
  error('Invalid PAIRIDX');
end
planeEqn = vertcat(ds.planeFit.plane);
planeEqn(:,4) = planeEqn(:,4)/job.metadata.pixelSize(1); % to pixels

% Coordinates.
if opts.trackIdx
  trackIdx = pairIdx;
else
  trackIdx = ds.sisterList(1).trackPairs(pairIdx,1:2);
end
featIdx = horzcat(ds.trackList(trackIdx).featIndx);
coords1 = nan(nFrames,3);
coords2 = nan(nFrames,3);
for i=1:nFrames
  if all(~isnan(featIdx(i,:)))
    coords1(i,:) = ds.initCoord(i).allCoordPix(featIdx(i,1),[1 2 3]);
    coords2(i,:) = ds.initCoord(i).allCoordPix(featIdx(i,2),[1 2 3]);
  end
end

% If 2D image set Z to 1.
if ~job.metadata.is3D
  coords1(:,3) = 1;
  coords2(:,3) = 1;
end

minXY = min([coords1(:,1:2); coords2(:,1:2)]);
maxXY = max([coords1(:,1:2); coords2(:,1:2)]);
maxDist = norm(maxXY-minXY,2);
width = (1+opts.margin)*maxDist; % Kymograph width in image pixels.
marginPix = opts.margin*maxDist;

% Kymograph RGB image array.
K = zeros(nFrames,opts.points,3);
irange = nan(3,2);
valid = zeros(nFrames,1);
plateDist = nan;
sisVec = [];
for i=1:nFrames
  if isnan(coords1(i,1)) && isempty(sisVec);
    continue;
  end
  valid(i) = 1;

  if ~isnan(coords1(i,1))
    % Sister-sister vector.
    s1 = coords1(i,1:2);
    s2 = coords2(i,1:2);
    sisVec = s1-s2;
    sisVec = sisVec / norm(sisVec,2);
  end

  % Endpoints.
  switch opts.fixture
    case 'sisters'
      % Sister centered.
      endPts = [s1 + marginPix*sisVec;
                s2 - marginPix*sisVec];
    case 'plate'
      % Fixed by plate.
      pn = planeEqn(i,1:2);
      pd = planeEqn(i,4);
      if isnan(plateDist)
        % Compute distance to plate from end of kymograph.
        plateDist = dot(pn,s1 + marginPix*sisVec) - pd;
      end

      % Maintain fixed distance to plane, find distance to move along sisVec.
      beta = (pd + plateDist - dot(s1,pn))/dot(sisVec,pn);
      endPts = [s1 + beta*sisVec; s1 - (width-beta)*sisVec];
  end

  for j=1:nChannels
    % Read frame. FIXME using sister 1 z plane arbitrarily
    plane = kitReadImagePlane(reader,md,i,j,coords1(i,3),job.crop);

    % First frame defines contrast stretch.
    if isnan(irange(mapChans(j),1))
      irange(mapChans(j),:) = stretchlim(plane,opts.saturate(j,:)/100);
    end

    K(i,:,mapChans(j)) = improfile(plane,endPts(:,1),endPts(:,2),...
                                   opts.points,'bilinear')';
  end
end

reader.close();

% Trim off untracked regions.
K = K(valid==1,:,:);

% Contrast stretch.
for j=1:3
  K(:,:,j) = imadjust(K(:,:,j), irange(j,:), []);
end

% If jet mode extract chosen channel and colormap.
if opts.jet > 0
  K = K(:,:,mapChans(opts.jet));
  n = 65536;
  K = gray2ind(K,n);
  K = ind2rgb(K,jet(n));
end

% Rescale.
K = imresize(K, 'Scale', fliplr(opts.scale), 'Method', opts.scaleMethod);

% Scale bar.
if ~isempty(opts.scaleBar)
  off = 5;
  wid = 4;
  len = opts.scale(1)*opts.scaleBar(1)*round(1/md.pixelSize(1));
  % um horizontal
  K(end-off-wid+1:end-off,off:off+len-1,:) = ones(wid,len,3);
  % sec vertical
  len = opts.scale(2)*opts.scaleBar(2);
  K(end-off-len+1:end-off,off:off+wid-1,:) = ones(len,wid,3);
end

% Convert NaNs to zero.
K(isnan(K)) = 0;

if opts.rotate
  K = imrotate(K,90);
end

figure;
clf;
imshow(K,[]);

% Save kymograph as PNG.
if ~isempty(opts.outfile)
  if isempty(strfind(opts.outfile,'.png'))
    opts.outfile = [opts.outfile '.png'];
  end
  imwrite(K,opts.outfile);
end
