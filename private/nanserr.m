function y = nanserr(x,dim)
% NANSERR Standard error of the mean, ignoring NaNs.
%   Y = NANSERR(X) returns the sample standard error of the mean of the values
%   in X, treating NaNs as missing values.  For a vector input, Y is the
%   standard deviation of the non-NaN elements of X.  For a matrix input, Y is a
%   row vector containing the standard error of the non-NaN elements in each
%   column of X. For N-D arrays, NANSERR operates along the first non-singleton
%   dimension of X.
%
%   Y = NANSERR(X,DIM) takes the standard deviation along dimension
%   DIM of X.
%
%   See also STD, NANVAR, NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.

sz = size(x);

if isempty(x), y = nan; return; end

if nargin < 2
  dim = find(sz ~= 1, 1);
  if isempty(dim), dim = 1; end
end

y = nanstd(x,0,dim);
n = sum(~isnan(x),dim);
y = y./sqrt(n);
