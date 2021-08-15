function spots=histcutSpotsSimplified(img,options,dataProperties,verbose)
% Find spots using histogram mode cutoff.
%
% Based on code from MaKi, by (mostly) K. Jaqaman and J. Dorn.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith and J. U. Harrison
% Copyright (c) 2019 J. U. Harrison

if nargin<4
  verbose = 0;
end

img = img/max(img(:));

[sx,sy,sz] = size(img);
if sz == 1
  ndims = 2;
else
  ndims = 3;
end

filters = createFilters(ndims,dataProperties);

% Filter image.
if verLessThan('images','9.2')
  imageF = fastGauss3D(img,filters.signalP(1:3),filters.signalP(4:6));
  background = fastGauss3D(img,filters.backgroundP(1:3),filters.backgroundP(4:6));
else
  imageF = imgaussfilt3(img,filters.signalP(1:3),'FilterSize',filters.signalP(4:6));
  background = imgaussfilt3(img,filters.backgroundP(1:3),'FilterSize',filters.backgroundP(4:6));
end

% get local maxima from the image
bw = imregionalmax(imageF);
localMax1DIndx = find(bw);
[x,y,z]=ind2sub([sx sy sz],localMax1DIndx);
locMax=[x,y,z];
if any(x>sx) || any(y>sy) || (ndims==3 && any(z>sz))
  error('invalid index');
end

% get signal strength of maxima
amp = imageF(localMax1DIndx) - background(localMax1DIndx);

% Filter out spots by cutting first histogram mode.
z = zeros(size(amp));
% amplitude
cutoff = splitModes(amp(~isApproxEqual(amp,z,[],'absolute')),[],[],[]);

if verbose
  figure(1)
  histogram(amp);
  hold on;
  yl = ylim;
  plot(repmat(cutoff,[1 2]),yl,'r');
  title('amplitude');
end

% Check if sufficient spots found in descending order of cutoff strictness.
nn = sum(amp>cutoff);
nn = (nn - options.minSpotsPerFrame)/(options.maxSpotsPerFrame-options.minSpotsPerFrame);
% If number of spots does not fall within range, then try a different threshold. Limit number of iterations of this
%
iter = 0;
realisticNumSpots = 0;
while (~realisticNumSpots && (iter < 5)) %TODO: set num iter as an option
if nn > 1
  %too many spots detected, increase threshold
  cutoff = cutoff*2;
  nn = sum(amp>cutoff);  
  nn = (nn - options.minSpotsPerFrame)/(options.maxSpotsPerFrame-options.minSpotsPerFrame);
elseif nn < 0
  %too few spots detected
  cutoff = cutoff/2;
  nn = sum(amp>cutoff);
  nn = (nn - options.minSpotsPerFrame)/(options.maxSpotsPerFrame-options.minSpotsPerFrame);
else
% Otherwise, cutoff gives acceptable number of spots, so use it
  realisticNumSpots=1;
end
iter = iter+1;
end

    passIdx = (amp > cutoff(1));

% keep these spots
spots = locMax(passIdx,:);

if verbose
  h = figure(2);
  imshow(max(img,[],3),[]);
  scale = 3;
  figpos = get(h,'Position');
  set(h,'Position',[figpos(1:2) figpos(3:4)*scale]);
  if ~isempty(spots)
    hold on;
    plot(spots(:,2),spots(:,1),'rx');
    hold off;
  end
  title('spots');
end
