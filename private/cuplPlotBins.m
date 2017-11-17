function cuplPlotBins(x,y,ye,varargin)
% CUPLPLOTBINS Draws bins by value
%
%   CUPLPLOTBINS(X,Y,YE) Draws Y binned into classes based on X, with errorbars
%   YE. Y and YE may be a vector for a single group of correlations or a matrix
%   for multiple correlations in columns.
%
%   CUPLPLOTBINS(X,Y,...) Options for plotting can be specified as
%   string/variable pairs, as follows:
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
plotYLabel = 'Value';

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

if ~all(size(y)==size(ye))
    error('Y and YE must have same size.');
end

% Make x column vector, if not.
if size(x,1)==1
    x=x';
end

hold on;
for i=1:size(y,2)
    errorbar(x,y(:,i),ye(:,i),['.' colours(i)]);
end
hold off;

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
    out = zeros(size(y,1),2*size(y,2));
    out(:,1:2:end) = y;
    out(:,2:2:end) = ye;
    dlmwrite(strrep(output,'pdf','csv'),[x,out],'precision',6);
end
