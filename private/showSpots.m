function showSpots(img,spots,pixels)
% Display spots on top of image.

if nargin<3
  pixels=0; % Set to rgb triple for colour.
end

h=figure(1);
if (isscalar(pixels) && pixels) || ~isscalar(pixels)
  % If 3D image, max project.
  img = max(img,[],3);

  % Make 3 layers out of original image (normalized).
  img = img/max(img(:));
  img = repmat(img,[1 1 3]);

  % Show pixels.
  for i=1:size(spots,1)
    img(spots(i,1),spots(i,2),:) = [0 0 1];
  end

  % Plot image.
  imshow(img);
  drawnow;
else
  imshow(max(img,[],3));
  hold on;
  if ~isempty(spots)
    plot(spots(:,2),spots(:,1),'rx');
  end
  hold off;
end

% Scale up to at least 512 pixels
figpos = getpixelposition(h);
if any(figpos(3:4)<256)
  set(h,'Units','Pixels','Position',[figpos(1:2) (512/min(figpos(3:4)))*figpos(3:4)]);
end
