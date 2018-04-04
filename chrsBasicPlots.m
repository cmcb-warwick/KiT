function chrsBasicPlots(zetaStructures,varargin)
% CHRSBASICPLOTS Produces plots of a given statistic of chromatic shift
% measurements.
%
%    CHRSBASICPLOTS(ZETASTRUCTURES,...) Outputs a figure representing the
%    distribution of chromatic shift measurements, collated in each cell
%    element of ZETASTRUCTURES, across multiple experiments.
%    ZETASTRUCTURES are produced by chrsZetaMeasurements.
%
%    Options, defaults in {}:-
%
%    filtered: 0 or {1}. Whether or not to use measurements filtered using
%       the parameters defined in job.options.chrShift.
%
%    legend: {'Expt 1', ...}, 'off', or similar. Names for each experiment.
%    
%    stat: {''}, or one of the following:
%           - 'zeta3D'
%           - 'zeta2D'
%           - 'zetaXYZ'
%       The statistic to be printed to screen. If no statistic is provided,
%       the user will be prompted.
%
% Copyright (c) 2017 C. A. Smith


% default options
opts.filtered = 1;
opts.legend = {''};
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

% process input
if ~iscell(zetaStructures)
  nExpts = 1;
  zetaStructures = {zetaStructures};
else
  nExpts = length(zetaStructures);
  % construct legend
  if isempty(opts.legend{1}) && nExpts>1
    for iExpt=1:nExpts
      opts.legend{iExpt} = ['Expt ' num2str(iExpt)];
    end
  end
end

% ask the user which stat they would like
statsList = {'zeta3D' ,'zeta2D' ,'zetaXYZ'};
if ~ismember(opts.stat,statsList)
    fprintf('Output statistics options:\n')
    for iStat=1:length(statsList)
        fprintf('    %i) %s\n',iStat,statsList{iStat});
    end
    prompt = sprintf('Please type the number statistic you would like: ');
    result = input(prompt);
    opts.stat = statsList{result};
end

switch opts.stat

    case 'zeta3D'
        
        if opts.filtered
          for iExpt = 1:nExpts
            zeta3D{iExpt} = zetaStructures{iExpt}.filtered.zeta.threeD(:)*1000;
          end
        else
          for iExpt = 1:nExpts
            zeta3D{iExpt} = zetaStructures.raw.zeta.threeD(:)*1000;
          end
        end
        
        if ~iscell(opts.legend) && strcmp(opts.legend,'off')
          title = sprintf('\\zeta_{3D}');
        else
          title = sprintf('\\zeta_{3D}: %s',opts.legend{1});
          for iExpt = 2:nExpts
            title = [title ' vs. ' opts.legend{iExpt}];
          end
        end
        compareHistograms(zeta3D,'axisLimits',[0 300],'title',title,'xLabel','\zeta_{3D} (nm)',...
            'legend',opts.legend);
        
    case 'zeta2D'
        
        if opts.filtered
          for iExpt = 1:nExpts
            zeta2D{iExpt} = zetaStructures{iExpt}.filtered.zeta.twoD(:)*1000;
          end
        else
          for iExpt = 1:nExpts
            zeta2D{iExpt} = zetaStructures.raw.zeta.twoD(:)*1000;
          end
        end
        
        if ~iscell(opts.legend) && strcmp(opts.legend,'off')
          title = sprintf('\\zeta_{2D}');
        else
          title = sprintf('\\zeta_{2D}: %s',opts.legend{1});
          for iExpt = 2:nExpts
            title = [title ' vs. ' opts.legend{iExpt}];
          end
        end
        compareHistograms(zeta2D,'axisLimits',[0 300],'title',title,'xLabel','\zeta_{2D} (nm)',...
            'legend',opts.legend);
        
    case 'zetaXYZ'
        
        if opts.filtered
          for iExpt = 1:nExpts
            zetaX{iExpt} = zetaStructures{iExpt}.filtered.zeta.x(:)*1000;
            zetaY{iExpt} = zetaStructures{iExpt}.filtered.zeta.y(:)*1000;
            zetaZ{iExpt} = zetaStructures{iExpt}.filtered.zeta.z(:)*1000;
          end
        else
          for iExpt = 1:nExpts  
            zetaX{iExpt} = zetaStructures{iExpt}.raw.zeta.x(:)*1000;
            zetaY{iExpt} = zetaStructures{iExpt}.raw.zeta.y(:)*1000;
            zetaZ{iExpt} = zetaStructures{iExpt}.raw.zeta.z(:)*1000;
          end
        end
        
        figure; clf
        if ~iscell(opts.legend) && strcmp(opts.legend,'off')
          title = sprintf('\\zeta_{x}');
        else
          title = sprintf('\\zeta_{x}: %s',opts.legend{1});
          for iExpt = 2:nExpts
            title = [title ' vs. ' opts.legend{iExpt}];
          end
        end
        subplot(1,3,1)
        compareHistograms(zetaX,'nBins',10,'axisLimits',[-300 300],...
            'title',title,'xLabel','\zeta_{x} (nm)','withinFig',1,...
            'legend',opts.legend);
        subplot(1,3,2)
        title(8) = 'y';
        compareHistograms(zetaY,'nBins',10,'axisLimits',[-300 300],...
            'title',title,'xLabel','\zeta_{y} (nm)','yLabel','','withinFig',1,...
            'legend',opts.legend);
        subplot(1,3,3)
        title(8) = 'z';
        compareHistograms(zetaZ,'nBins',10,'axisLimits',[-300 300],...
            'title',title,'xLabel','\zeta_{z} (nm)','yLabel','','withinFig',1,...
            'legend',opts.legend);
        
    otherwise
        
        error('Statistic requested is not yet built into this version of chrsBasicPlots. See later release.')
    
end

end