function m = cuplCellSerr(analysis, x)
% CUPLCELLSERR Compute serr averages per cell for quantity x
%
%   M = CUPLCELLSERR(ANALYSIS,X) Compute standard error of the mean M per cell
%   for quantity X
%
%   X NxM quantity to compute standard error for. Computation is done across
%   columns, result is NxP, where P is number of cells. NaNs are ignored via
%   nanserr.
%
% Copyright (c) 2013 Jonathan Armond

% Allocate storage.
n = analysis.nCells;
m = zeros(size(x,1),n);

for i=1:n
  % Indexes of accepted sisters in accepted cell i.
  sisterIdx = analysis.sisterCellIdx == i;
  % Compute cell serr.
  m(:,i) = nanserr(x(:,sisterIdx),2);
end
