function compareHistograms(data,varargin)
%COMPAREHISTOGRAMS Produces overlaid histograms of up to 5 distributions.
%
%    COMPAREHISTOGRAMS(DATA,...) For n distributions of data in cell form
%    {nx1}, histograms are plotted with adaptable level of transparency.
%    Number of bins, and range over which data is plotted can be adjusted.
%    Fits can be made to the histogram distributions, peak values given.
%
%    Example of use:-
%    
%    compareHistograms({inputData1,inputData2},'nBins',50,'fitType','Normal')
%        will plot a single distribution of data with 50 bins, and attempt
%        to fit a normal distribution.
%
%    Options, defaults in {}:-
%
%    axisLimits: {[]} or two numbers. The range over which the data is
%             plotted. The second of the two numbers must be larger than
%             the first.
%
%    transparency: {0.5} or number between 0 and 1. The transparency of the
%             histograms, so that transparency increases the closer the
%             number is to 1.
%
%    fitDistribution: {'none'} or one of the following:
%             - 'Beta'
%             - 'Gamma'
%             - 'Normal'
%             - 'tLocationScale'
%             The normalised distribution which is best fit to the data, a
%             peak value given. Normally used to calculate most common
%             value that the data holds. NB: distributions are
%             case-sensitive, i.e. 'gamma' or 'tLocationScale ' will be
%             rejected.
%
%    legend: {'Expt 1','Expt 2','Expt 3','Expt 4','Expt 5'} or similar.
%             Titles for each experiment, as will be shown in the legend,
%             and best fit peak values.
%    
%    nBins: {20} or number. Number of bins over which to plot the data.
%
%    plotStyle: 'bars' or {'outline'}. Will either plot filled bars, or the
%             outline of the bars.
%
%    title: {''} or a string. Title of the plot.
%
%    withinFig: {0} or 1. Defines whether or not to plot the histogram
%             within the current figure environment.
%
%    xLabel: {''} or a string. Label for the x-axis.
%
%    yLabel: {'normalised frequency'} or similar. Label for the y-axis.
%   
% Created by: 2013 Christopher Smith

% Set defaults
opts.axisLimits = [];
opts.transparency = 0.5;
opts.fitDistribution = 'none'; % can also use 'Beta','Gamma','Normal','tLocationScale'
opts.legends = {'Expt 1','Expt 2','Expt 3','Expt 4', 'Expt 5'};
opts.nBins = 20;
opts.plotStyle = 'outline'; % can also use 'bars'
opts.title = '';
opts.withinFig = 0;
opts.xLabel = '';
opts.yLabel = 'normalised frequency';

% Process options
opts = processOptions(opts, varargin{:});
CC = [ 0 , 0 , 1;
       1 , 0 , 0;
      0.5, 1 , 0;
      0.5, 0 , 1;
       0 , 0 , 0];
nCols = size(CC,1);
       

%% Data handling

% Process data
if iscell(data)
    nExpts = length(data);
    for iExpt = 1:nExpts
        % Remove any unwanted data
        data{iExpt} = data{iExpt}(~isnan(data{iExpt}));
        if ~isempty(opts.axisLimits)
            data{iExpt} = data{iExpt}(data{iExpt}<opts.axisLimits(2));
            data{iExpt} = data{iExpt}(data{iExpt}>opts.axisLimits(1));
        end
        emptyData(iExpt) = isempty(data{iExpt});
        % Ensure all datasets are column vectors
        data{iExpt} = data{iExpt}(:);
    end
    if sum(emptyData) == nExpts
        error('Empty dataset provided.')
    end
else
    % Ensure data is in a cell in column vector form (if it isn't empty)
    % Remove any NaNs in the process
    if isempty(data)
        error('Empty dataset provided.')
    else
        % Remove any unwanted data
        data = data(~isnan(data));
        if ~isempty(opts.axisLimits)
            data = data(data<opts.axisLimits(2));
            data = data(data>opts.axisLimits(1));
        end
        % Ensure all datasets are column vectors
        data = data(:);
        data = {data};
        nExpts=1; emptyData=0;
    end
end

% Check have enough colours for the number of experiments
if nExpts > nCols
   CC = repmat(CC,ceil(nExpts/nCols),1); 
end

% check legends
if ~iscell(opts.legends) && strcmp(opts.legends,'off')
  legOff = 1;
elseif iscell(opts.legends)
  while length(opts.legends)<nExpts
      legLength = length(opts.legends);
      opts.legends{end+1} = ['Expt ' legLength+1];
  end
  legOff = 0;
else
  warning('Legend provided not in correct format. Converting to default.')
  opts.legends = {'Expt 1','Expt 2','Expt 3','Expt 4', 'Expt 5'};
  legOff = 0;
end

% Find raw minimum and maximum of full dataset
minVal = nanmin(cat(1,data{:}));
maxVal = nanmax(cat(1,data{:}));
% Find whether positive or negative values dominate the range of values,
% and define that as the normalising factor
if abs(minVal) > abs(maxVal)
    normFact = abs(minVal);
else
    normFact = abs(maxVal);
end
% Normalise minimum and maximum values
normMinVal = minVal/normFact;
normMaxVal = maxVal/normFact;
% Calculate the range, and therefore bin separations and bin centres
range = ceil(abs(normMaxVal-normMinVal));
normBinsep = range/opts.nBins;
normBins = floor(normMinVal)+(normBinsep/2) : normBinsep : ceil(normMaxVal)-(normBinsep/2);
% De-normalise bins
bins = normBins*normFact;

for iExpt = 1:nExpts
    % Calculate histogram bin intensities.
    hExpt{iExpt} = hist(data{iExpt},bins);
    hExpt{iExpt} = hExpt{iExpt}/trapz(bins,hExpt{iExpt});
    
    % Get number of datapoints and adjust legends.
    nPoints = length(data{iExpt});
    opts.legends{iExpt} = sprintf('%s (n=%i)',opts.legends{iExpt},nPoints);
end

%% Plotting data

if ~opts.withinFig
  figure()
  clf
end
hold on

switch opts.plotStyle
    
    case 'bars' % plots histograms as individual bars
        
        for iExpt = 1:nExpts
            bExpt{iExpt} = bar(bins,hExpt{iExpt},1,'FaceColor',CC(iExpt,:));
        end
        
    case 'outline'
        
        % duplicate bins to collate corners of boxes
        corners = [bins-normBinsep/2*normFact, bins(end)+normBinsep/2*normFact];
        corners = repmat(corners,2,1);
        corners = corners(:)';
        
        % duplicate data accordingly and buffer each end with zeros
        for iExpt=1:nExpts
            if ~isempty(hExpt{iExpt})
                hOutline{iExpt} = hExpt{iExpt};
                hOutline{iExpt} = repmat(hOutline{iExpt},2,1);
                hOutline{iExpt} = hOutline{iExpt}(:)';
                hOutline{iExpt} = [0, hOutline{iExpt}, 0];
                plot(corners,hOutline{iExpt},'Color',CC(iExpt,:),'LineWidth',2)
            end 
        end
        
end

%% Aesthetics

faceAlpha = 1-opts.transparency;
% Adjust bar colours and transparencies
if strcmp(opts.plotStyle,'bars')
    for iExpt = 1:nExpts
        set(bExpt{iExpt}, 'faceAlpha', faceAlpha, 'FaceColor', CC(iExpt,:))
    end
end

% Adjust x-axis limits
if isempty(opts.axisLimits)
    opts.axisLimits = [min(bins) max(bins)];
end
xlim(opts.axisLimits);

% Apply legend, title and axis labels
set(gca,'FontSize',20)

if ~legOff
  legend(opts.legends{~emptyData}); legend boxoff
end
title(opts.title, 'FontSize', 20)
ylabel(opts.yLabel)
xlabel(opts.xLabel)

end


