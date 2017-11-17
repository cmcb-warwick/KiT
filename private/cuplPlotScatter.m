function cuplPlotScatter(x,y,varargin)
%CUPLPLOTCORRELATION  Draws scatter plot of correlations by distance
%
%   CUPLPLOTCORRELATION(X,Y) Draws scatter plot of correlations Y at distance
%   X. X and Y may be a vector for a single group of correlations or a cell
%   array of vectors for multiple groups of correlations.
%
%   CUPLPLOTCORRELATION(X,Y,...) Options for plotting can be specified
%   as string/variable pairs, as follows:
%
%       'Colours': Sequence of line colours (String; default 'krbgmcy').
%
%       'PlotTitle': Add title to top of graph (String).
%
%       'PlotYlim': Limit y-axis (2-element float vector).
%
%       'Output': Save plot to file (String).
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

% Defaults.
colours = 'krbgmcy';
plotTitle = '';
plotYLim = [];
output = [];
plotXLabel = 'Distance (\mum)';
plotYLabel = 'Correlation';

% Read options.
options = {'colours','plotTitle','plotYLim','output','plotXLabel','plotYLabel'};
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

if ~iscell(y)
    y={y};
end
if ~iscell(x)
    x={x};
end
if size(y)~=size(x)
    error('X and Y must have same number of cell elements.');
end

hold on;
for i=1:length(y)
    plot(x{i},y{i},['.' colours(i)]);
end

if ~isempty(plotYLim)
    ylim(plotYLim);
end

% Draw axis line.
axisLine = line(xlim(gca),[0 0]);
set(axisLine,'LineStyle',':');
axisAnnotation = get(axisLine,'Annotation');
axisLegendEntry = get(axisAnnotation,'LegendInformation');
set(axisLegendEntry,'IconDisplayStyle','off');

% Label graph.
xlabel(plotXLabel);
ylabel(plotYLabel);
title(plotTitle);

% Annotation.
%text(0.05*max(xlim),-0.85,['n = ' num2str(size(nonsis_dist,1))],'FontSize',8);

% Save.
if ~isempty(output)
    saveas(gcf,output);
    for i=1:length(y)
        if i==1
            dlmwrite(strrep(output,'pdf','csv'),[x{i},y{i}],'precision',6);
        else
            dlmwrite(strrep(output,'pdf','csv'),[x{i},y{i}],'-append',...
                     'precision',6,'roffset',2);
        end
    end
end
