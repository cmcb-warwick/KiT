function kitPlayMovie(filename,outfile,docrop)
% KITPLAYMOVIE Play movie file

if nargin<1 || isempty(filename)
  [filename,pathname] = uigetfile(...
    kitSupportedFormats(1),'Locate movie to play');
  if isempty(filename)
    disp('No file selected');
    return
  end
else
  [pathname,filename,ext] = fileparts(filename);
  filename = [filename ext];
end
filename = fullfile(pathname,filename);

if ~exist(filename,'file')
  error('File not found');
end

if nargin<2
  outfile=[];
end

if nargin>2 && docrop == 1
  [crop,cropSize] = kitCropMovie(filename);
  if mod(cropSize(1),2) ~= 0
    crop(3) = crop(3)-1;
    cropSize(1) = cropSize(1)-1;
  end
  if mod(cropSize(2),2) ~= 0
    crop(4) = crop(4)-1;
    cropSize(2) = cropSize(2)-1;
  end
else
  crop = [];
  cropSize = [];
end

% Open movie.
[md,reader] = kitOpenMovie(filename);
opts.saturate = [1 99.9];
if size(opts.saturate,1)<md.nChannels
  opts.saturate = repmat(opts.saturate,[md.nChannels,1]);
end
mapChans = [2 1 3]; % Green, red, blue
maxMergeChannels = 3;
dt = md.frameTime(1,2)-md.frameTime(1,1);
len = md.frameTime(1,end)-md.frameTime(1,1);
fprintf('dt = %.2fs  length = %.1f min\n', dt, len/60);

figure;
clf;
if ~isempty(outfile)
  vWriter = VideoWriter(outfile, 'MPEG-4');
  vWriter.FrameRate = 5;
  vWriter.Quality = 95;
  open(vWriter);
end

if isempty(cropSize)
  imgSize = md.frameSize(1:2)
else
  imgSize = cropSize(1:2)
end

for i=1:md.nFrames
  rgbImg = zeros([imgSize, 3]);
  for c=1:min(md.nChannels, maxMergeChannels)
    % Read stack.
    img = kitReadImageStack(reader, md, i, c, crop, 0);

    % Max project.
    img = max(img, [], 3);

    % First frame defines contrast stretch.
    if i==1
      irange(c,:)=stretchlim(img,opts.saturate(c,:)/100);
    end

    % Contrast stretch.
    rgbImg(:,:,mapChans(c)) = imadjust(img, irange(c,:), []);
  end

  imshow(rgbImg);
  axis image;
  drawnow;
  if ~isempty(outfile)
    % Save frame.
    writeVideo(vWriter, getframe);
  end
end

if ~isempty(outfile)
  close(vWriter);
end
