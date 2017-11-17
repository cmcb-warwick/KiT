function cuplPlotCorrelation(t,y,ye,varargin)
%CUPLPLOTCORRELATION  Draws correlation plot
%
%   CUPLPLOTCORRELATION(T,Y,YE) Draws plot of correlations Y at times T with
%   confidence interval YE. Y and YE may be a vector for a single correlation or
%   matrix for multiple.
%
%   CUPLPLOTCORRELATION(T,Y,YE,...) Options for plotting can be specified
%   as string/variable pairs, as follows:
%
%       'Colours': Sequence of line colours (String; default 'krbgmcy').
%
%       'PlotTitle': Add title to top of graph (String).
%
%       'PlotYlim': Limit y-axis (2-element float vector).
%
%       'MaxTime': Limit x-axis to [0,MaxTime] (Float).
%
%       'Output': Save plot to file (String).
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

% Defaults.
colours = 'krbgmcy';
maxTime = 0;
plotTitle = '';
plotYLim = [];
output = [];

% Read options.
options = {'colours','maxTime','plotTitle','plotYLim','output'};
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

% Make t column vector, if not.
if size(t,1)==1
    t=t';
end

if ~isempty(ye) && ~all(size(y)==size(ye))
        error('Y and YE must be same size');
end

hold on;
for i=1:size(y,2)
    plot(t,y(:,i),['-' colours(i)]);
    if ~isempty(ye)
        plot(t,y(:,i)+ye(:,i),[':' colours(i)]);
        plot(t,y(:,i)-ye(:,i),[':' colours(i)]);
    end
end

if maxTime>0
    xlim([0 maxTime]);
else
    xlim([0 max(t)]);
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
hold off;

% Label graph.
xlabel('Time lag (s)');
ylabel('Correlation');
title(plotTitle);

% Save.
if ~isempty(output)
    saveas(gcf,output);
    if isempty(ye)
        out=[t,y];
    else
        out=[t,y,ye];
    end
    dlmwrite(strrep(output,'pdf','csv'),out,'precision',6);
end
