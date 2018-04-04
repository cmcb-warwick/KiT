function [figH,rgbImg,z]=kitMovieProj(movieFileName,zPlane,roi,progress)
% KITMOVIEPROJ Presents movie projection in Z and T
%
%  SYNOPSIS kitShowMovieProj(movieFileName)
%
%  INPUT movieFileName: Filename of movie to crop.
%        zPlane: Show single z-plane (optional).
%
% Copyright (c) 2012 Jonathan W. Armond

if nargin<2
  zPlane=[];
end
if nargin<3
  roi=[];
end
if nargin<4
  progress=0;
end

% Open movie.
if progress
  h = waitbar(0,'Opening movie');
end
[metadata,reader] = kitOpenMovie(movieFileName);

% Load subset of frames and overlay.
nImages = 10;
frameList = unique(round(1:metadata.nFrames/nImages:metadata.nFrames));

% Loop to create max projections.
maxMergeChannels = 3;
imgSz = [metadata.frameSize(1:2), 3];

centreSz = 0.5; % Take central percentage for locating threshold.
border = (1-centreSz)/2;
cx(1) = round(imgSz(1)*border);
cx(2) = round(cx(1) + imgSz(1)*centreSz);
cy(1) = round(imgSz(2)*border);
cy(2) = round(cy(1) + imgSz(2)*centreSz);

rgbImg = zeros(imgSz);
mapChan = [2 1 3];
for c=1:min([maxMergeChannels, metadata.nChannels])
  maxProj = zeros(metadata.frameSize(1:2));
  for f=1:length(frameList)
    if progress
      waitbar(f/length(frameList),h);;
    end
    % Read stack.
    if isempty(zPlane)
      img = kitReadImageStack(reader, metadata, frameList(f), c);
    else
      img = kitReadImagePlane(reader, metadata, frameList(f), c, zPlane);
    end
    z = size(img,3);
    % Average max projection of each frame.
    maxProj = maxProj + max(img,[],3);
  end
  % Normalize.
  maxProj = maxProj/length(frameList);

  % Merge into RGB image.
  centreImg = maxProj(cx(1):cx(2),cy(1):cy(2),:);
  lb = splitModes(centreImg(centreImg>0));
  if isempty(lb)
    % Can't identify background. Use sensible default.
    lb = 0.75;
  end
  slim = stretchlim(centreImg,[lb 1]);
  rgbImg(:,:,mapChan(c)) = imadjust(maxProj,slim,[]);
end

if progress
  close(h);
end
figH=figure;
imshow(rgbImg);

if ~isempty(roi)
  rectangle('Position',roi,'LineWidth',2,'EdgeColor','y');
end

