function heatMap(A,B,varargin)

if nargin<2 || isempty(B) || isempty(A)
    error('Please supply two vectors of data for comparison.')
elseif ~isvector(A) || ~isvector(B)
    error('Input data must be in vector form.')
elseif length(A) ~= length(B)
    error('Two vectors of equal length required for comparison.')
end

% set default options
opts.colormap = 'jet';
opts.enforceMax = [];
opts.invData = 0;
opts.nBins = 100;
opts.title = '';
opts.withinFig = 0;
opts.xLabel = '';
opts.xLimits = [];
opts.yLabel = '';
opts.yLimits = [];

% process options
opts = processOptions(opts,varargin{:});

% enforce input to be column vectors
A = A(:); B = B(:);
% combine the data
data = [A,B];

% process x- and y-limits if necessary
xLimits = opts.xLimits; yLimits = opts.yLimits;
if isempty(xLimits)
    xLimits = [floor(min(A)) ceil(max(A))];
end
if isempty(yLimits)
    yLimits = [floor(min(B)) ceil(max(B))];
end

% take only the data within these limits
data = data(data(:,1)>xLimits(1),:);
data = data(data(:,1)<xLimits(2),:);

data = data(data(:,2)>yLimits(1),:);
data = data(data(:,2)<yLimits(2),:);

% find bin size based on limits, then define space required for hist3
binSizeX = (xLimits(2) - xLimits(1))./opts.nBins;
binSizeY = (yLimits(2) - yLimits(1))./opts.nBins;

histSpace{1} = xLimits(1):binSizeX:xLimits(2);
histSpace{2} = yLimits(1):binSizeY:yLimits(2);

% perform bivariate histogram analysis
h = hist3(data,histSpace);
h = h';
h = h/trapz(h(:));

if ~isempty(opts.enforceMax)
    cLim = [0 opts.enforceMax];
else
    cLim = [min(h(:)) max(h(:))];
end

if opts.invData
    yLimits = fliplr(yLimits);
end

% produce image
if ~opts.withinFig
    figure()
end
image(xLimits,yLimits,h,'CDataMapping','scaled')
colorbar
set(gca,'YDir','normal','FontSize',20,'CLim',cLim)
colormap(opts.colormap)
xlabel(opts.xLabel,'FontSize',20); ylabel(opts.yLabel,'FontSize',20)
title(opts.title,'FontSize',20);

end