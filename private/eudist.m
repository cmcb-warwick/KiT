function d=eudist(x,y)
% EUDIST Calculate Euclidean distance between two points
%
%   D = EUDIST(X,Y) Calculate Euclidean distance between two points. Ignores
%   NaNs via nansum.
%
%   X, Y can be Nx3 matrices, in which case the distance between the points
%   specified by each corresponding row is returned as a column vector D length
%   N.

% If vectors, ensure row.
if isvector(x) && size(x,1) > 1
  x = x';
end
if isvector(y) && size(y,1) > 1
  y = y';
end

if size(x,2) ~= 3 || size(y,2) ~= 3
  error(['X and Y must be either both vectors of length 3, or one Nx3 ' ...
         'matrix and a vector of length 3, or both Nx3 matrices']);
end

d = sqrt(sum(bsxfun(@minus,x,y).^2, 2));
