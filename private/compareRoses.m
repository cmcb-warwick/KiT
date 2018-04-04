function f = compareRoses(data,varargin)
%COMPAREROSES Produces overlaid rose plots of up to 3 distributions.
%
%    COMPAREROSES(DATA,...) For n distributions of data in cell form {nx1},
%    roses are plotted with adaptable level of transparency. Number of bins
%    can be adjusted.
%
%    Example of use:-
%    
%    compareRoses({inputData1,inputData2},'nBins',24,'plotStyle','outline')
%        will plot a single distribution of data with 24 bins, with only a
%        line plotted around the edge of the boxed data.
%
%    Options, defaults in {}:-
%
%    transparency: {0.5} or number in [0 1]. The transparency of the
%             histograms, so that transparency increases the closer the
%             number is to 1.
%
%    legend: {'Expt 1','Expt 2','Expt 3','Expt 4','Expt 5'} or similar.
%             Titles for each experiment, as will be shown in the legend.
%    
%    nBins: {12} or number. Number of bins over which to plot the data.
%
%    plotStyle: {'bars'}, 'outline' or 'border'. Will either plot filled bars, or the
%             outline of the bars.
%
%    title: {''} or a string. Title of the plot.
%
%    scaling: {'normalise'} or 'maximum'. Will either normalise to set area
%             below to 1, or 'maximum' will scale so that the maximum bi
%             is 1.
%
%    units: {'deg'} or 'rad'. Units of the input data.
%
%    withinFig: {0} or 1. Defines whether or not to plot the histogram
%             within the current figure environment.
%
%    xLabel: {''} or a string. Label for the x-axis.
%   
% Created by: 2015 Christopher Smith

%% Pre-processing

% Set defaults
opts.transparency = 0.5;
opts.legend = {'Expt 1','Expt 2','Expt 3'};
opts.nBins = 12;
opts.plotStyle = 'outline'; % can also be 'bars' or 'border'
opts.title = '';
opts.scaling = 'maximum'; % can also be 'normalise'
opts.units = 'deg'; % can also be 'rad'
opts.withinFig = 0;
opts.xLabel = '';

% Process options
opts = processOptions(opts, varargin{:});
CC = [ 1 , 0 , 0;
       0 ,0.5, 0;
      0.5, 1 , 0;
      0.5, 0 , 1;
       0 , 0 , 0];
nCols = size(CC,1);

% Process data
if iscell(data)
    nExpts = length(data);
    for iExpt = 1:nExpts
        data{iExpt} = data{iExpt}(~isnan(data{iExpt}));
        emptyData(iExpt) = isempty(data{iExpt});
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
        data = data(~isnan(data));
        data = data(:);
        data = {data};
        nExpts=1; emptyData=0;
    end
end

% Check have enough colours for the number of experiments
if nExpts > nCols
   CC = repmat(CC,ceil(nExpts/nCols),1); 
end

% Ensure data is in radians for rose plotting
if strcmp(opts.units,'deg')
    for iExpt = 1:nExpts
        data{iExpt} = data{iExpt}*(pi/180);
    end
end

% Find bins
binSep = 2*pi/opts.nBins;
bins = 2*pi/opts.nBins:binSep:2*pi;

% Calculate histogram bin intensities
for iExpt = 1:nExpts
    if ~isempty(data{iExpt})
        [the{iExpt} rho{iExpt}] = rose(data{iExpt},bins);
    end
end

% Calculate distribution based on normalisation factor
switch opts.scaling
    case 'normalise'
        for iExpt = 1:nExpts
            if ~isempty(data{iExpt})
                rho{iExpt} = rho{iExpt}./trapz(rho{iExpt});
            end
        end
    case 'maximum'
        for iExpt = 1:nExpts
            if ~isempty(data{iExpt})
                rho{iExpt} = rho{iExpt}./max(rho{iExpt});
            end
        end
end

%% Plotting data

if ~opts.withinFig
    figure();
end

% Pre-plot rose to get plotting environment, then delete the plot
for iExpt = 1:nExpts
    if ~isempty(data{iExpt})
        h{iExpt} = polar(the{iExpt},rho{iExpt});
        the{iExpt} = get(h{iExpt},'XData');
        rho{iExpt} = get(h{iExpt},'YData');
        delete(h{iExpt})
    end
end
hold on

switch opts.plotStyle
    
    case 'bars' % plots histograms as individual bars
        
        faceAlpha = 1 - opts.transparency;
        % plot first bins first to ensure legend is correct
        for iExpt = 1:nExpts
            hp = patch('XData',the{iExpt}(1:4),'YData',rho{iExpt}(1:4));
            set(hp,'FaceColor',CC(iExpt,:),'FaceAlpha',faceAlpha,'LineStyle','none');
        end
        
        % plot remaining bins
        for iExpt = 1:nExpts
            if ~isempty(data{iExpt})
                for iBin = 2:opts.nBins
                    hp = patch('XData',the{iExpt}(4*iBin-3:4*iBin),'YData',rho{iExpt}(4*iBin-3:4*iBin));
                    set(hp,'FaceColor',CC(iExpt,:),'FaceAlpha',faceAlpha,'LineStyle','none');
                end
                plot(the{iExpt},rho{iExpt},'k','LineWidth',2)
            end
        end
        
    case 'border'
        
        for iExpt = 1:nExpts
            if ~isempty(data{iExpt})
                
                nNodes = length(the{iExpt})/4;
                
                the{iExpt} = reshape(the{iExpt},4,nNodes)';
                the{iExpt} = nanmean(the{iExpt}(:,2:3),2);
                the{iExpt}(end+1) = the{iExpt}(1);
                
                rho{iExpt} = reshape(rho{iExpt},4,nNodes)';
                rho{iExpt} = nanmean(rho{iExpt}(:,2:3),2);
                rho{iExpt}(end+1) = rho{iExpt}(1);
                
                plot(the{iExpt},rho{iExpt},'Color',CC(iExpt,:),'LineWidth',2)
                
            end
        end
        
    case 'outline'
        
        for iExpt = 1:nExpts
            if ~isempty(data{iExpt})
                
                nNodes = length(the{iExpt})/4;
                
                the{iExpt} = reshape(the{iExpt},4,nNodes)';
                the{iExpt} = the{iExpt}(:,2:3);
                the{iExpt} = reshape(the{iExpt}',1,nNodes*2);
                the{iExpt}(end+1) = the{iExpt}(1);
                
                rho{iExpt} = reshape(rho{iExpt},4,nNodes)';
                rho{iExpt} = rho{iExpt}(:,2:3);
                rho{iExpt} = reshape(rho{iExpt}',1,nNodes*2);
                rho{iExpt}(end+1) = rho{iExpt}(1);
                
                plot(the{iExpt},rho{iExpt},'Color',CC(iExpt,:),'LineWidth',2)
            end
        end
        
end

%% Aesthetics
set(gca,'FontSize',20)
xlabel(opts.xLabel)
title(opts.title)
legend(opts.legend{~emptyData}); legend boxoff

end