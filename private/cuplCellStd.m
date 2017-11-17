function m = cuplCellStd(analysis, x, flag)
% CUPLCELLSTD Compute standard deviation per cell for quantity x
%
%   M = CUPLCELLSTD(ANALYSIS,X,FLAG) Compute standard deviation M per cell for
%   quantity X.
%
%   X NxM quantity to compute standard deviation for. Computation is done
%    across columns, result is NxP, where P is number of cells. NaNs are ignored
%    via nanstd.
%
%   FLAG Optional, is passed through to nanstd and has the same meaning.
%
% Copyright (c) 2013 Jonathan Armond

if nargin<3
  flag = 0;
end

% Allocate storage.
n = analysis.nCells;
m = zeros(size(x,1),n);

for i=1:n
  % Indexes of accepted sisters in accepted cell i.
  sisterIdx = analysis.sisterCellIdx == i;
  % Compute cell std.
  m(:,i) = nanstd(x(:,sisterIdx),flag,2);
end
