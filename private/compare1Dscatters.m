function compare1Dscatters(data,varargin)
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
%    legend: {'Expt 1', ...} or similar. Titles for each experiment, as
%             will be shown on the x-axis beneath each distribution.
%    
%    title: {''} or a string. Title of the plot.
%
%    withinFig: {0} or 1. Defines whether or not to plot the 1D scatter
%             within the current figure environment.
%
%    yLabel: {'measurement (units)'} or similar. Label for the y-axis.
%   
% Copyright (c) 2017 C.A.Smith

% Set defaults
opts.axisLimits = [];
opts.legend = [];
opts.title = '';
opts.withinFig = 0;
opts.yLabel = 'measurement (units)';

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
    
% Ensure that there are legends for each experiment
if isempty(opts.legend)
  for iExpt = 1:nExpts
    legends{iExpt} = ['Expt ' num2str(iExpt)];
  end
else
  legends = opts.legend;
end
  

% Find raw minimum and maximum of full dataset
minVal = nanmin(cat(1,data{:}));
maxVal = nanmax(cat(1,data{:}));

%% Plotting data

if ~opts.withinFig
    figure()
    clf
end
hold on

for iExpt = 1:nExpts
  
  nPoints = length(data{iExpt});
  legends{iExpt} = sprintf('%s (n=%i)',legends{iExpt},nPoints);
  
  xData = ones(nPoints,1)*iExpt - (rand(nPoints,1)-0.5)/5;
    
  scatter(xData,data{iExpt},'o','MarkerEdgeColor',CC(iExpt,:))
  
end

%% Aesthetics

% Adjust y-axis limits
if ~isempty(opts.axisLimits)
  ylim(opts.axisLimits)  
end
xlim([0 nExpts+1])
    
% Apply legend, title and axis labels
set(gca,'FontSize',20)

title(opts.title, 'FontSize', 20)
ylabel(opts.yLabel)
set(gca,'XTick',1:nExpts)
set(gca,'XTickLabel',legends)

end
