function kitMakeSpotMovie(job, varargin)
% KITMAKESPOTMOVIE Renders a movie with marked spots
%
%    KITMAKESPOTMOVIE(JOB,...) Renders the movie associated with JOB with
%    various annotations. Supply options as string/value pairs following JOB.
%
%    Options, defaults in {}:-
%
%    outfile: {none} or filename. Filename to write movie to.
%
%    rotate: {0} or 1. Transforms image to align with fitted coordinate system.
%
%    plotPlane: {0} or 1. Plot axes of fitted coordinate system (X - red, Y - cyan).
%
%    plotTracks: {0} or 1. Draw spot tracks.
%
%    plotSpots: 0 or {1}. Draw spots.
%
%    plotSisterFrames: {[]} or cell array. Plot arrows next to sisters in
%    specific frames. Cell array with sister indexs in first element and
%    corresponding frame lists in subsequent entrys.
%
%    annotate: {0} or 1. Add annotation with frame time, number etc..
%
%    showAttach: {0} or 1. Indicate pole attachment with triangle.
%
%    slow: {0} or number of seconds. Delay between displaying each frame to
%    allow easier viewing.
%
%    trackChannel: {first tracked channel} or channel number. Tracked channel
%    to get spots from.
%
%    zoomTrack: {0} or track number. Zoom into a particular track(s).
%
%    zoomMargin: {0.1} or number. Enlarge track boundary region by factor
%    when using zoomTrack.
%
%    intensityGraph: {0} or 1. Draw animated intensity graph for zoomed
%    track.
%
%    centerZoom: {0} or 1. When zoomTrack is enabled, if centerZoom is 1,
%    then the zoomed spot will be central, otherwise the movie will be zoomed
%    to the extent of the spots motion.
%
%    plotMaxCoord: {0} or channel number. Plot location of maximum intensity
%    pixel within intensity mask, when zoomTrack enabled.
%
%    drawMask: {0} or 1. Draw approximation of intensity mask, when zoomTrack
%    enabled.
%
%    codec: {'Motion JPEG AVI'} Codec to use for image compression. See help
%    VideoWriter.
%
%    saturate: {[1 100]} or 2 element vector or scalar, between 0 and 100. As a
%    scalar, percentage of pixel values, at either end of histogram, to
%    saturate. Increase to enhance constrast more. Alternatively, supply a
%    vector [LOW HIGH] to specify percentage to saturate at low and high pixels
%    values. NB 2 is equivalent to [1 99]. Or supply a matrix to adjust
%    independently for each channel.
%
%    channelMap: {[2 1 3]} or a 3-element vector containing integers between
%    1 and 3. Controls which channels are shown in which colors. 1 is red, 2
%    is green and 3 is blue.
%
%    normalize: {0} or 1. Normalize image values to range [0,1]. Useful for some
%    floating-point data with unknown range, since MATLAB expects [0,1].
%
%    transpose: {0} or 1. Transpose image data.
%
%    scale: {1} or number. Scale image.
%
%    See also: KITPLAYMOVIE, KITMAKEKYMOGRAPH.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<1
  error('Must supply JOB');
end

% Set defaults.
opts.rotate = 0;
opts.plotPlane = 0;
opts.zoomTrack = 0;
opts.zoomMargin = 0.1;
opts.intensityGraph = 0;
opts.centerZoom = 0;
opts.plotTracks = 0;
opts.plotMaxCoord = 0; % Channel to plot max coord for.
opts.plotSisterFrames = [];
opts.drawMask = 0;
opts.annotate = 1; % Add annotation text.
opts.plotSpots = 1;
opts.showAttach = 0; % Show attached pole with triangle.
opts.trackChannel = 0; % Default to first tracked channel.
opts.codec = 'Motion JPEG AVI';
opts.outfile = [];
opts.slow = 0;
opts.saturate = [1 100];
opts.channelMap = [2 1 3]; % Green, red, blue
opts.normalize = 0;
opts.transpose = 0;
opts.fliplr = 0; % Flip x and y spot coordinates.
opts.labelSpots = 0;
opts.plotSisters = 0;
opts.scale = 1;

% Process options.
opts = processOptions(opts, varargin{:});

% Handle options.
if opts.plotSisters == 1
  opts.plotSpots = 0;
end

colors = presetColors();

% Open movie.
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie));

h = figure;
clf;
if ~isempty(opts.outfile)
  vWriter = VideoWriter(opts.outfile, opts.codec);
  vWriter.FrameRate = 5;
  vWriter.Quality = 95;
  open(vWriter);
end

coordSysChan = job.options.coordSystemChannel;
if opts.trackChannel == 0
  % Plot first non-empty tracking channel
  plotChans = find(~cellfun(@isempty,job.dataStruct),1);
else
  plotChans = opts.trackChannel;
end
initCoord = job.dataStruct{plotChans}.initCoord;
tracks = job.dataStruct{plotChans}.tracks;
sisterList = job.dataStruct{plotChans}.sisterList;
trackList = job.dataStruct{plotChans}.trackList;
if isfield(job.dataStruct{plotChans},'trackInt')
  trackInt = job.dataStruct{plotChans}.trackInt;
end
nTracks = length(trackList);
if opts.fliplr
  for i=1:nTracks
    trackList(i).coords(:,1:2) = trackList(i).coords(:,[2 1]);
  end
end

% Extract pixel track coordinates.
for i=1:nTracks
  trackList(i).pixCoords = nan(md.nFrames,6);
  for t=1:md.nFrames
    featIndx = trackList(i).featIndx(t);
    if ~isnan(featIndx)
      trackList(i).pixCoords(t,:) = ...
          job.dataStruct{plotChans}.initCoord(t).allCoordPix(featIndx,:);
    end
  end
  if opts.fliplr
    trackList(i).pixCoords(:,1:2) = trackList(i).pixCoords(:,[2 1]);
  end
end

% Inset graph parameters.
insetSz = 0.5;
insetOff = [1 -0.1];

% marker size.
markerSize = ceil(job.options.maskRadius / job.metadata.pixelSize(1));

if opts.zoomTrack > 0
  track = trackList(opts.zoomTrack);
  minCoords = [inf inf];
  maxCoords = [-inf -inf];
  for i=1:length(opts.zoomTrack)
    % Find max extent of zoom track.
    coords = nan(md.nFrames,2);
    for j=1:md.nFrames
      if ~isnan(track(i).coords(j,1))
        cds = track(i).coords(j,1:2);
        cds = cds ./ job.metadata.pixelSize(1:2);
        cds = cds + job.cropSize([2 1])/2;
        coords(j,:) = cds;
      end
    end
    minCoords = nanmin([minCoords; coords]);
    maxCoords = nanmax([maxCoords; coords]);
  end

  maxDim = max(maxCoords-minCoords)*(1+opts.zoomMargin);
  centreCoords = (maxCoords+minCoords)/2;
  extent = [centreCoords-maxDim centreCoords+maxDim]; % bottom-left x,y top-right x,y
  centreFrame = 20;
  % Assign direction of first track.
  trackDir = [track(1).direction; 0];
end


% Read frame by frame.
mapChans = opts.channelMap;
%markers = 'xo+*sd';
maxMergeChannels = 3;
if opts.transpose
  dims = [2 1];
else
  dims = [1 2];
end

if size(opts.saturate,1)<md.nChannels
  opts.saturate = repmat(opts.saturate,[md.nChannels,1]);
end

for i=1:md.nFrames
  rgbImg = zeros([job.cropSize(dims), 3]);
  for c=1:min(md.nChannels, maxMergeChannels)
    % Read stack.
    img = kitReadImageStack(reader, md, i, c, job.crop, opts.normalize);

    % Max project.
    img = max(img, [], 3);
    if opts.transpose
      img = img';
    end

    % First frame defines contrast stretch.
    if i==1
      irange(c,:)=stretchlim(img,opts.saturate(c,:)/100);
    end

    % Contrast stretch.
    rgbImg(:,:,mapChans(c)) = imadjust(img, irange(c,:), []);
  end

  sz = size(rgbImg);
  imgCentre = sz([2 1])/2;
  % Transform image.
  planeFit = job.dataStruct{coordSysChan}.planeFit(i);
  origin = planeFit.planeOrigin(1:2)./job.metadata.pixelSize(1:2);
  if opts.rotate && ~isempty(planeFit.planeVectors)
    if isfield(planeFit, 'tform')
      rgbImg = imtransform(rgbImg, planeFit.tform,...
                           'XData',[1 sz(2)],'YData',[1 sz(1)],'Size',[sz(1) sz(2)]);
    else
      % Translate to image coordinates (0,0).
      transMat = eye(3);
      transMat(3,1:2) = -origin;
      % Rotate into alignment with coordinate system.
      rotMat = eye(3);
      rotMat(1:2,1:2) = planeFit.planeVectors(1:2,1:2);
      % Translate to image centre.
      invTransMat = eye(3);
      invTransMat(3,1:2) = imgCentre;
      % Composite transform => translate, rotate, translate back.
      tform = maketform(...
        'composite', maketform('affine', invTransMat), ...
        maketform('affine', rotMat), maketform('affine', transMat));
      rgbImg = imtransform(...
        rgbImg,tform,'XData',[1 sz(2)],'YData',[1 sz(1)],...
        'Size',[sz(1) sz(2)]);
    end
  end
  imshow(rgbImg);
  axis image;

  hold on
  if opts.plotTracks == 1
    % Plot tracks.
    for c=plotChans
      if opts.zoomTrack == 0
        range = 1:nTracks;
      else
        range = opts.zoomTrack;
      end
      for k=range
        % Extract corresponding coordinates.
        coords = job.dataStruct{c}.trackList(k).coords(1:i,1:2);
        coords = coords ./ repmat(job.metadata.pixelSize(1:2),i,1);
        % Rows <=> Y, Cols <=> X.
        coords = coords + repmat(imgCentre,i,1);
        plot(coords(:,1), coords(:,2), 'w-');
      end
    end
  elseif opts.plotSpots == 1
    % Plot spots.
    for c=plotChans
      if opts.rotate
        trackCoords = horzcat(trackList.coords);
        coords = [trackCoords(i,1:6:end)',trackCoords(i,2:6:end)'];
        coords = coords ./ repmat(job.metadata.pixelSize(1:2),size(coords,1),1);
        % Rows <=> Y, Cols <=> X.
        coords = coords + repmat(imgCentre,size(coords,1),1);
      else
        trackCoords = horzcat(trackList.pixCoords);
        coords = [trackCoords(i,1:6:end)',trackCoords(i,2:6:end)'];
      end

      if ~isempty(coords)
        hold on;
        if opts.zoomTrack == 0
          range = 1:size(coords,1);
        else
          range = opts.zoomTrack;
        end

        % Sort in pole attachments.
        attach = vertcat(trackList(range).attach);
        if opts.showAttach
          % Decided by sister.
          plot(coords(range(attach==-2),1), coords(range(attach==-2),2), 'w<');
          plot(coords(range(attach==+2),1), coords(range(attach==+2),2), 'w>');
          % Decided by location.
          plot(coords(range(attach==-1),1), coords(range(attach==-1),2), 'g<');
          plot(coords(range(attach==+1),1), coords(range(attach==+1),2), 'g>');
        else
          if opts.labelSpots
            for k=1:length(range)
              text(coords(range(k),1),coords(range(k),2),num2str(range(k)), ...
                   'color',[1 1 1]);
            end
          else
            plot(coords(range,1),coords(range,2),'wx');
          end
        end
      end
    end
  elseif opts.plotSisters == 1 ||  ~isempty(opts.plotSisterFrames)
    % Plot sisters.
    for c=plotChans
      sisterCoords1 = horzcat(sisterList.coords1);
      sisterCoords2 = horzcat(sisterList.coords2);
      coords1 = [sisterCoords1(i,1:6:end)',sisterCoords1(i,2:6:end)'];
      coords2 = [sisterCoords2(i,1:6:end)',sisterCoords2(i,2:6:end)'];
      coords1 = coords1 ./ repmat(job.metadata.pixelSize(1:2),size(coords1,1),1);
      coords2 = coords2 ./ repmat(job.metadata.pixelSize(1:2),size(coords2,1),1);
      % Rows <=> Y, Cols <=> X.
      coords1 = coords1 + repmat(imgCentre,size(coords1,1),1);
      coords2 = coords2 + repmat(imgCentre,size(coords2,1),1);

      if opts.plotSisters
        for k=1:size(coords1,1)
          plot(coords1(k,1),coords1(k,2),'x','Color',colors(k,:))
          plot(coords2(k,1),coords2(k,2),'x','Color',colors(k,:))
        end
      else
        pairs = sisterList(1).trackPairs(:,1:2);
        idx = opts.plotSisterFrames{1};
        frames = opts.plotSisterFrames(2:end);
        for k=1:length(idx)
          if ismember(i,frames{k})
            % get pixel coords.
            p = pairs(idx(k),1:2);
            if length(tracks(p(1)).tracksFeatIndxCG)>=i
              f1 = tracks(p(1)).tracksFeatIndxCG(i);
              if f1>0
                plot(initCoord(i).allCoordPix(f1,1),initCoord(i).allCoordPix(f1,2),'rx')
              end
            end
            if length(tracks(p(2)).tracksFeatIndxCG)>=i
              f2 = tracks(p(2)).tracksFeatIndxCG(i);
              if f2>0
                plot(initCoord(i).allCoordPix(f2,1),initCoord(i).allCoordPix(f2,2),'rx')
              end
            end
            %plot(coords1(idx(k),1),coords1(idx(k),2),'rx')
            %plot(coords2(idx(k),1),coords2(idx(k),2),'rx')
          end
        end
      end
    end
  end

  if opts.drawMask
    for c=plotChans
      % Represent pixel mask.
      rotMat = eye(3);
      rotMat(1:2,1:2) = planeFit.planeVectors(1:2,1:2);
      for j=1:length(trackInt)
        % Transform mask coords.
        maskCoord = trackInt(j).maskCoord(i,1:2);
        if isnan(maskCoord(1))
          continue;
        end
        maskCoord = [(maskCoord - origin) 1]*rotMat;
        maskCoord = imgCentre + maskCoord(1:2);

        if strcmp(job.options.maskShape,'semicircle')
          drawSemicircle(maskCoord(1),maskCoord(2),markerSize,trackList(j).attach,'w');
        else
          drawCircle(maskCoord(1),maskCoord(2),markerSize,'w');
        end
      end
    end
  end

  c=coordSysChan;
  % Plot plane fit.
  if opts.plotPlane && ~isempty(job.dataStruct{c}.planeFit(i).planeVectors)
    px = [job.dataStruct{c}.planeFit(i).planeVectors(1:2,1)' 0];
    py = [job.dataStruct{c}.planeFit(i).planeVectors(1:2,2)' 0];
    po = [job.dataStruct{c}.planeFit(i).planeOrigin(1:2)./job.metadata.pixelSize(1:2) 1];
    if opts.rotate
      po = [imgCentre 1];
      % Plotting occurs in normal coordinate system.
      rotMat = eye(3);
      rotMat(1:2,1:2) = planeFit.planeVectors(1:2,1:2);
      px = px * rotMat;
      py = py * rotMat;
    end
    plength = 40;
    p1 = [po; po+plength*px];
    p2 = [po; po+plength*py];
    plot(p1(:,1), p1(:,2), '-r'); % X
    plot(p2(:,1), p2(:,2), '-c'); % Y
                                  %plot(po(1), po(2),'yo');
  end

  % Plot max coordinate.
  if opts.plotMaxCoord > 0
    c = opts.plotMaxCoord;
    rotMat = eye(3);
    rotMat(1:2,1:2) = planeFit.planeVectors(1:2,1:2);
    for j=1:length(trackInt)
      % Transform max coords.
      maxCoord = trackInt(j).maxCoord(i,3*(c-1)+1:3*c-1);
      if isnan(maxCoord(1))
        continue;
      end
      maxCoord = [(maxCoord - origin) 1]*rotMat;
      maxCoord = round(imgCentre + maxCoord(1:2));
      plot(maxCoord(1),maxCoord(2),'m+');
    end
  end


  % Zoom plot
  if opts.zoomTrack > 0
    if opts.centerZoom == 0
      xlim(extent([1 3]));
      ylim(extent([2 4]));
    else
      % Centre frame on spot.
      if ~isnan(coords(i,1))
        cy = floor(coords(i,1));
        cx = floor(coords(i,2));
      end
      ylim([cx-centreFrame cx+centreFrame]);
      xlim([cy-centreFrame cy+centreFrame]);
    end
  end

  % Plot intensity graph.
  if (isvector(opts.zoomTrack) || opts.zoomTrack > 0) && opts.intensityGraph > 0
    c = opts.intensityGraph;
    j = opts.zoomTrack;
    xl = xlim();
    yl = ylim();
    graphSz = insetSz*[xl(2)-xl(1), yl(2)-yl(1)];
    graphOff = insetOff.*graphSz + [xl(1),yl(1)];
    [gx,gy] = transformPoints(job.metadata.frameTime,trackInt(j).intensity(:,c),...
                              graphSz,graphOff);
    plot(gx(1:i),yl(2)-yl(1)+gy(1:i),'w-');
  end

  % Draw timestamp.
  if opts.annotate == 1
    t = job.metadata.frameTime(1,i); % Take first z-slice as timestamp.
    sec = floor(t);
    msec = floor(1000*(t-sec));
    tstamp = sprintf('T %03d.%03d F %03d',sec,msec,i);
    if opts.zoomTrack > 0
      % Add P/AP direction.
      if trackDir(i)>0
        dir = 'P';
      elseif trackDir(i)<0
        dir = 'AP';
      else
        dir = 'N';
      end
      tstamp = [tstamp ' D ' dir];
      if opts.plotMaxCoord > 0
        c = opts.plotMaxCoord;
        % Add angle/dist to max.
        zTrack = trackInt(opts.zoomTrack);
        tstamp = [tstamp sprintf(' MAX % +4.0f,% 3.0f',...
                                 zTrack.angleToMax(i,c)*180/pi,zTrack.distToMax(i,c)*1000)];
      end
    end
    text(1,6,tstamp,'units','pixels','color',[1 1 1]);
  end


  hold off

  if opts.scale ~= 1 && i == 1
    figpos = get(h,'Position');
    set(h,'Position',[figpos(1:2) figpos(3:4)*opts.scale]);
  end

  if ~isempty(opts.outfile)
    % Save frame.
    writeVideo(vWriter, getframe);
  end

  if opts.slow > 0
    pause(opts.slow);
  else
    drawnow;
  end
end

  if ~isempty(opts.outfile)
    close(vWriter);
  end
end

%% LOCAL FUNCTIONS

function drawCircle(x,y,r,color)
% Draws circle.

% Estimate pixels in circumference.
c = 2*pi*r;
theta = linspace(0,2*pi,ceil(c));
cx = x + r*cos(theta);
cy = y + r*sin(theta);
plot(cx, cy, [color '-']);

end

function drawSemicircle(x,y,r,dir,color)
% Draws y-axis aligned semi-circle.

% Estimate pixels in circumference.
c = 2*pi*r;
theta = linspace(pi/2,3*pi/2,ceil(c));
if dir>0
  theta = theta + pi;
end
cx = x + r*cos(theta);
cy = y + r*sin(theta);
cx(end+1) = cx(1);
cy(end+1) = cy(1);
plot(cx, cy, [color '-']);

end

function [x,y]=transformPoints(x,y,sz,off)
% Transforms points for plotting into image coordinate system.

xlims = [nanmin(x),nanmax(x)];
ylims = [nanmin(y),nanmax(y)];
x = ((x-xlims(1))/xlims(2)) * sz(1) + off(1);
y = ylims(2) - ((y-ylims(1))/ylims(2)) * sz(2) + off(2);

end

function colors=presetColors()
  colors = [
         0         0    1.0000;
    1.0000         0         0;
         0    1.0000         0;
         0         0    0.1724;
    1.0000    0.1034    0.7241;
    1.0000    0.8276         0;
         0    0.3448         0;
    0.5172    0.5172    1.0000;
    0.6207    0.3103    0.2759;
         0    1.0000    0.7586;
         0    0.5172    0.5862;
         0         0    0.4828;
    0.5862    0.8276    0.3103;
    0.9655    0.6207    0.8621;
    0.8276    0.0690    1.0000;
    0.4828    0.1034    0.4138;
    0.9655    0.0690    0.3793;
    1.0000    0.7586    0.5172;
    0.1379    0.1379    0.0345;
    0.5517    0.6552    0.4828;
    0.9655    0.5172    0.0345;
    0.5172    0.4483         0;
    0.4483    0.9655    1.0000;
    0.6207    0.7586    1.0000;
    0.4483    0.3793    0.4828;
    0.6207         0         0;
         0    0.3103    1.0000;
         0    0.2759    0.5862;
    0.8276    1.0000         0;
    0.7241    0.3103    0.8276;
    0.2414         0    0.1034;
    0.9310    1.0000    0.6897;
    1.0000    0.4828    0.3793;
    0.2759    1.0000    0.4828;
    0.0690    0.6552    0.3793;
    0.8276    0.6552    0.6552;
    0.8276    0.3103    0.5172;
    0.4138         0    0.7586;
    0.1724    0.3793    0.2759;
         0    0.5862    0.9655];
end
