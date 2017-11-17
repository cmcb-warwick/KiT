function [bins,edges,centres,m,sd,se]=binData(x,y,n,range)
%BINDATA  Bin data Y by X and return bin means.
%
%   BINS = BINDATA(X,Y) Bins data Y by X into 10 bins and returns the 
%   datapoints in each bin in a cell array.
%
%   BINDATA(X,Y,N) Specify the number of bins in N.
%
%   BINDATA(X,Y,N,RANGE) Specify lower and upper range of bins in vector RANGE
%   in the form [LOWER,UPPER].
%
%   [BINS,EDGES,CENTRES,M,SD,SE] = BINDATA(X,Y) Also returns the right bin
%   EDGES, the bin CENTRES, and the MEAN, standard deviations (SD) and standard
%   errors (SE) of the datapoints.
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

% Error checking.
if ~isvector(x) || ~isvector(y)
    error('X and Y must be vectors.');
end

% Defaults.
if nargin<3
    n = 10;
end

if nargin<4
    range = [min(x) max(x)];
end

% Allocate storage.
m = zeros(n,1);
if nargout >= 4
    sd = zeros(n,1);
end
if nargout >= 5
    se = zeros(n,1);
end

% Bin edges.
edges = linspace(range(1),range(2),n+1);

% Bin X values.
[h, idx] = histc(x,edges);

bins = cell(n,1);
for i=1:n
    % Find datapoints in each bin.
    bins{i} = y(idx==i);
    
    % Take means, etc.
    if nargout >= 4
        m(i) = nanmean(bins{i});
    end
    if nargout >= 5
        sd(i) = nanstd(bins{i});
    end
    if nargout >= 6
        se(i) = sd(i)/sqrt(sum(~isnan(bins{i})));
    end
end

if nargout >= 3
    % Bin centres.
    centres = edges(1:end-1)' + (edges(2)-edges(1))/2;
end
