function a=autocorrelation(x,maxLags)
% AUTOCORRELATION  Calculate autocorrelation for timeseries x.
%
%   A = AUTOCORRELATION(X,MAXLAGS) Calculates autocorrelation A for timeseries
%   x, at lags of [-MAXLAGS,MAXLAGS].
%
%   X is a timeseries as a vector, or multiple timeseries as a matrix of columns.
%
%   MAXLAGS Optional, specify how many lags to calculate. Defaults to N/2,
%           where N is number of rows of X.
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

if size(x,1) == 1
  % x is row vector.
  x = x';
end

[n,m] = size(x);
if nargin < 2
  maxLags = floor(n/2);
end

% Means.
mx = nanmean(x);

% Autocorrelate each column.
a = zeros(maxLags+1,m);

% Calculate each column.
for lag=0:maxLags
  vec = (x(1:n-lag,:)-repmat(mx,n-lag,1)) .* ...
         (x(1+lag:n,:)-repmat(mx,n-lag,1));
  a(lag+1,:) = nanmean(vec);
end

% Normalize.
a = a./repmat(a(1,:),maxLags+1,1);
