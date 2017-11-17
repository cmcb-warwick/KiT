function cuplPlotHistogram(x,varargin)
%CUPLPLOTHISTOGRAM  Draws correlation plot
%
%   CUPLPLOTHISTOGRAM(X) Draws histogram of values X.
%
%   CUPLPLOTHISTOGRAM(X,...) Options for plotting can be specified
%   as string/variable pairs, as follows:
%
%       'PlotTitle': Add title to top of graph (String).
%
%       'PlotYlim': Limit y-axis (2-element float vector).
%
%       'NumBins': Number of bins for histogram (Integer; Default 10).
%
%       'Output': Save plot to file (String).
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

% Defaults.
plotTitle = '';
plotYLim = [];
output = [];
numBins = [];
plotXLabel = 'Distance (\mum)';
plotYLabel = 'Correlation';

% Read options.
options = {'plotTitle','plotYLim','output','numBins','plotXLabel','plotYLabel'};
for i=1:2:length(varargin)
    name = lower(varargin{i});
    value = varargin{i+1};
    matches = strcmp(lower(options),name);
    if all(matches==0)
        error(['Unknown option: ' name]);
    end
    idx = find(matches,1);
    eval([options{idx} '=value;']);
end

% Create figure.
h = figure;
if isempty(numBins)
  hg = histogram(x);
else
  hg = histogram(x,numBins);
end

if ~isempty(plotYLim)
    ylim(plotYLim);
end

% Label graph.
xlabel(plotXLabel);
ylabel(plotYLabel);
title(plotTitle);

% Annotate graph.
fontSize = 12;
th = 0.05;
tx = 0.05;
ty = 0.95;
text(tx,ty,sprintf('mean = %.4f',nanmean(x)),'FontSize',fontSize,'Units','Normalized');
ty=ty-th;
text(tx,ty,sprintf('median = %.4f',nanmedian(x)),'FontSize',fontSize,'Units','Normalized');
ty=ty-th;
text(tx,ty,sprintf('std = %.4f',nanstd(x)),'FontSize',fontSize,'Units','Normalized');
ty=ty-th;
text(tx,ty,sprintf('n = %d',sum(~isnan(x))),'FontSize',fontSize,'Units','Normalized');


% Save.
if ~isempty(output)
    saveas(gcf,output);
    dlmwrite(strrep(output,'pdf','csv'),[xout; n]','precision',6);
end
