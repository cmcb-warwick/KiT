function cols=indexToCols(idx,dim)
% INDEXTOCOLS Return columns indexes for given track/sister indexes.
%
%   COLS = INDEXTOCOLS(IDX,DIM) Return columns indexes for given track/sister
%   indexes.
%
%   IDX Sister/track indexes.
%
%   DIM Optional. Specify 1, 2, or 3, or a vector of same, for x, y, z to return
%   only the indexes of those columns.
%
%   e.g. indexToCols([1 3]) => [1 2 3 7 8 9]
%        indexToCols([1 3],2) => [2 8]
%        indexToCols([1 3],[1 2]) => [1 2 7 8]

idx = idx-1;

if nargin<2
  dim = [1 2 3];
end

cols = [];
for i = 1:length(dim)
  cols = [cols, idx*3+dim(i)];
end

cols = sort(cols);
