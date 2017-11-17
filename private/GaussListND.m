function gaussList = GaussListND(coordList,sigma,center,intNorm)
% GAUSSLISTND calculates the value of a N-D Gaussian at specific pixel/voxel coordinates
%
% SYNOPSIS gaussList = GaussListND(coordList,sigma,center,intNorm)
%
% INPUT    coordList : m-by-n list of coordinates, where m is the number of
%                      coordinates and n the number of dimensions
%          sigma     : 1-by-n (or scalar): sigma of Gaussian
%          center    : (opt) 1-by-n vector of center of Gaussian.
%                      Default: zeros(1,n)
%          intNorm   : (opt) switch for how the Gaussian should be normed
%                      Default: 0
%                      0 - no norming. Max of Gaussian == 1
%                      1 - normed so that integral of infinite Gaussian = 1
%
% OUTPUT   gaussList : m-by-1 list of intensities. Intensity is the
%                      integral of the Gaussian over the pixel/voxel
%
% REMARKS  The code assumes that a pixel has the edge length 1!
%
% Copyright (c) 2005 Jonas Dorn

% check dimensionality of coordList. 
if isempty(coordList)
  error('you have to supply a list of coordinates for GaussList23D')
else
  [nCoords,nDims] = size(coordList);
end

% sigma
if ~isscalar(sigma) && length(sigma) ~= nDims;
  error('sigma has to be a scalar or a 1-by-n vector!')
end

% center
if nargin < 3 || isempty(center)
  center = zeros(nCoords,nDims);
else
  if ~isscalar(center) && length(center) ~= nDims;
    error('center has to be a scalar or a 1-by-n vector!')
  end
end

% intNorm
if nargin < 4 || isempty(intNorm)
  intNorm = 0;
end

%======================
% CALC GAUSSLIST
%======================

% 0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the
% pixel at 1 of a Gaussian with mean 0 and sigma 1

% convert coordList to 0/1
% (coordList - center)./sigma
coordList = bsxfun(@rdivide,bsxfun(@minus,coordList,center),sigma);

% double coordList as preparation for erfc
% cat(3,coordList-0.5./sigma, coordList+0.5./sigma)
coordList = cat(3,bsxfun(@minus,coordList,0.5./sigma),...
                bsxfun(@plus,coordList,0.5./sigma));

% calculate gaussList
gaussList = diff(0.5 * erfc(-coordList/sqrt(2)),1,3);
gaussList = prod(gaussList,2);

% norm gaussList
switch intNorm
  case 0
    if isscalar(sigma)
      sigmaProd = sigma^3;
    else
      sigmaProd = prod(sigma);
    end
    % "un-norm" Gaussian
    gaussList = gaussList*((2*pi)^(0.5*nDims)*sigmaProd);
  case 1
    % gaussList is already normed
end
