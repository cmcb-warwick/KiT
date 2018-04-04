function compareBars(data,varargin)
%COMPAREBARS Produces bar plots of up to 5 distributions.
%
%    COMPAREBARS(DATA,...) For n distributions of DATA in cell form
%    {nx1}, bars are plotted with error bars. Bars can either represent the
%    mean or median, while error bars can represent standard deviation,
%    standard error of the mean, or be omitted altogether.
%
%    Example of use:-
%    
%    compareBars({inputData1,inputData2},'stat','mean','errorBars','stdErr')
%        will plot two bars representing the mean of the data, and error
%        bars representing the standard error.
%
%    Options, defaults in {}:-
%
%    errorBars: {'stdDev'}, 'stdErr' or 'none'. Which statistic, if any, to
%             use to represent error bars.
%
%    legend: {'Expt 1', ...} or similar. Titles for each experiment, as
%             will be shown on the x-axis beneath each distribution.
%
%    stat: {'median'} or 'mean'. Which statistic to plot.
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
opts.errorBars = 'stdDev';
opts.legend = [];
opts.stat = 'median';
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
        % Ensure all datasets are column vectors
        data{iExpt} = data{iExpt}(:);
        emptyData(iExpt) = isempty(data{iExpt});
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
  

%% Plotting data

if ~opts.withinFig
    figure()
    clf
end
hold on

for iExpt = 1:nExpts
  
  xData = iExpt;
  nPoints = length(data{iExpt});
  legends{iExpt} = sprintf('%s (n=%i)',legends{iExpt},nPoints);
  
  switch opts.stat
    case 'mean'
      yValue = nanmean(data{iExpt});
    case 'median'
      yValue = nanmedian(data{iExpt});
    otherwise
      error('Statistic not recognised: %s.',opts.stat)
  end
  
  switch opts.errorBars
    case 'stdDev'
      errorBar = nanstd(data{iExpt});
    case 'stdErr'
      errorBar = nanserr(data{iExpt});
    case 'none'
      errorBar = 0;
    otherwise
      errorBar = 0;
      warning('Error statistic not recognised: %s.',opts.errorBars)
  end
    
  bar(xData,yValue,'FaceColor',CC(iExpt,:))
  if errorBar~=0
    line([xData xData],[yValue-errorBar yValue+errorBar],'Color','k')
    line([xData-0.1 xData+0.1],[yValue-errorBar yValue-errorBar],'Color','k')
    line([xData-0.1 xData+0.1],[yValue+errorBar yValue+errorBar],'Color','k')
  end
  
end

%% Aesthetics

% Adjust y-axis limits
xlim([0 nExpts+1])
    
% Apply legend, title and axis labels
set(gca,'FontSize',20)

title(opts.title, 'FontSize', 20)
ylabel(opts.yLabel)
set(gca,'XTick',1:nExpts)
set(gca,'XTickLabel',legends)

end
