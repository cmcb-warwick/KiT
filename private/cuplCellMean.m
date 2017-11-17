function m = cuplCellMean(analysis, x)
% CUPLCELLMEAN Compute mean averages per cell for quantity x
%
%   M = CUPLCELLMEAN(ANALYSIS,X) Compute mean averages M per cell for quantity X
%
%   X NxM quantity to average. Averaging is done across columns, result is
%   NxP, where P is number of cells. NaNs are ignored via nanmean.
%
% Copyright (c) 2013 Jonathan Armond

% Allocate storage.
n = analysis.nCells;
m = zeros(size(x,1),n);

for i=1:n
  % Indexes of sisters in cell i.
  sisterIdx = analysis.sisterCellIdx == i;
  % Compute cell mean.
  m(:,i) = nanmean(x(:,sisterIdx),2);
end
