function [F,J] = fitNGaussians3D(x0,img,index,psfSigma)
% FITNGAUSSIANS3D yields F, the difference between an image and a theoretical
% image produced by N Gaussians, and J, the Jacobian of F.
%
% SYNOPSIS [F,J] = FITNGAUSSIANS3D(X0,IMG,INDEX,PSFSIGMA)
%
% INPUT  X0      : initial guess of PSF positions and amplitudes and
%                  background noise.
%        IMG     : Image part being analyzed.
%        INDEX   : x,y,z-indices of pixels considered.
%        PSFSIGMA: Standard deviation of point spread function (in pixels).
%
% OUTPUT F       : Residuals from fitting an image with supplied
%                  Gaussians.
%        J       : The Jacobian matrix of F.
%
% REMARKS F = model image - real image, important to know if the sign of the
% residuals matters.
%
% Copyright (c) 2005 K. Jaqaman
% Copyright (c) 2012 Ed Harry
% Copyright (c) 2013 Jonathan W. Armond

% Extract background intensity from x0 and remove from vector.
bgAmp = x0(end);
x0(end) = [];

% Get number of PSFs considered.
numPSF = length(x0)/4;

% Reshape 4nx1 vector x0 into nx4 matrix.
x0 = reshape(x0,4,numPSF)';

% Extract PSF center positions and amplitudes.
psfPos = x0(:,1:3);
psfAmp = x0(:,4);

% Find minimum and maximum pixel indices.
minIndxX = min(index(:,1));
maxIndxX = max(index(:,1));
minIndxY = min(index(:,2));
maxIndxY = max(index(:,2));
minIndxZ = min(index(:,3));
maxIndxZ = max(index(:,3));

% Determine the contribution of each PSF (assuming amplitude 1) to a
% pixel based on its x-coordinate (needed to calculate F & J)
psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF);
for i=1:numPSF
    psfIntegX(:,i) = GaussListND((minIndxX:maxIndxX)',...
        psfSigma(1),psfPos(i,1));
end

% Determine the contribution of each PSF (assuming amplitude 1) to a
% pixel based on its y-coordinate (needed to calculate F & J).
psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF);
for i=1:numPSF
    psfIntegY(:,i) = GaussListND((minIndxY:maxIndxY)',...
        psfSigma(1),psfPos(i,2));
end

% Determine the contribution of each PSF (assuming amplitude 1) to a
% pixel based on its z-coordinate (needed to calculate F & J).
psfIntegZ = zeros(maxIndxZ-minIndxZ+1,numPSF);
for i=1:numPSF
    psfIntegZ(:,i) = GaussListND((minIndxZ:maxIndxZ)',...
        psfSigma(2),psfPos(i,3));
end

% Get xy-indices relative to minimum.
relIndxX = index(:,1) - minIndxX + 1;
relIndxY = index(:,2) - minIndxY + 1;
relIndxZ = index(:,3) - minIndxZ + 1;

% Calculate the value of F at all pixels.
F = bsxfun(@plus,...
           (sum(bsxfun(@times,psfAmp, psfIntegX(relIndxX,:)' .* ...
                       psfIntegY(relIndxY,:)' .* psfIntegZ(relIndxZ,:)'),1))',...
           bgAmp) - img;


if nargout > 1
  h = psfSigma(1)/1000;

  % Calculate the value of each PSF (assuming amplitude 1) at the
  % x-coordinates of the corners of all pixels (needed to calculate J).
  psfValueX = zeros(2*(maxIndxX-minIndxX+1),numPSF);
  for i=1:numPSF
    psfValueX(:,i) = GaussListND([(minIndxX:maxIndxX)-h (minIndxX:maxIndxX)+h]',psfSigma(1),psfPos(i,1));
  end

  % Calculate the value of each PSF (assuming amplitude 1) at the
  % y-coordinates of the corners of all pixels (needed to calculate J).
  psfValueY = zeros(2*(maxIndxY-minIndxY+1),numPSF);
  for i=1:numPSF
    psfValueY(:,i) = GaussListND([(minIndxY:maxIndxY)-h (minIndxY:maxIndxY)+h]',psfSigma(1),psfPos(i,2));
  end

  % Calculate the value of each PSF (assuming amplitude 1) at the
  % z-coordinates of the corners of all pixels (needed to calculate J).
  psfValueZ = zeros(2*(maxIndxZ-minIndxZ+1),numPSF);
  for i=1:numPSF
    psfValueZ(:,i) = GaussListND([(minIndxZ:maxIndxZ)-h (minIndxZ:maxIndxZ)+h]',psfSigma(2),psfPos(i,3));
  end

  % Get number of pixels in image.
  numPixel = length(img);

  % Calculate the derivative at all pixels.
  J = ones(numPixel,4*numPSF+1); % (last column for background amplitude)
  J(:,1:4:4*numPSF) = bsxfun(@times, psfAmp', ...
      (psfValueX(relIndxX,:) - psfValueX(maxIndxX-minIndxX+relIndxX+1,:))/(2*h) .* ...
      psfIntegY(relIndxY,:) .* psfIntegZ(relIndxZ,:)); %w.r.t. x
  J(:,2:4:4*numPSF) = bsxfun(@times, psfAmp', ...
      (psfValueY(relIndxY,:) - psfValueY(maxIndxY-minIndxY+relIndxY+1,:))/(2*h) .* ...
      psfIntegX(relIndxX,:) .* psfIntegZ(relIndxZ,:)); %w.r.t. y
  J(:,3:4:4*numPSF) = bsxfun(@times, psfAmp', ...
      (psfValueZ(relIndxZ,:) - psfValueZ(maxIndxZ-minIndxZ+relIndxZ+1,:))/(2*h) .* ...
      psfIntegX(relIndxX,:) .* psfIntegY(relIndxY,:)); %w.r.t. z
  J(:,4:4:4*numPSF) = psfIntegX(relIndxX,:) .* ...
      psfIntegY(relIndxY,:) .* psfIntegZ(relIndxZ,:); %w.r.t. amp
end
