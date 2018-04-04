function chrsBasicStats(zetaStructure,varargin)
% CHRSBASICSTATS Gives a list of statistics for delta measurements.
%
%    CHRSBASICSTATS(COMPILEDZETA,...) Calculates the median, mean and
%    standard error and deviation of measurements of zeta collated in
%    ZETASTRUCTURE and outputs to the command line. ZETASTRUCTURE is
%    produced by chrsZetaMeasurements.
%
%    Options, defaults in {}:-
%
%    filtered: 0 or {1}. Whether or not to give intensity-, neighbour- and
%       region-filtered measurements of zeta measurements.
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
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

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
          zeta3D = zetaStructure.filtered.zeta.threeD(:);
        else
          zeta3D = zetaStructure.raw.zeta.threeD(:);
        end
        
        median = nanmedian(zeta3D)*1000;
        mean   = nanmean(zeta3D)*1000;
        stdErr = nanserr(zeta3D)*1000;
        stdDev = nanstd(zeta3D)*1000;
        n      = min(sum(~isnan(zeta3D)));

        fprintf('\n3D zeta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'zeta2D'
        
        if opts.filtered
          zeta2D = zetaStructure.filtered.zeta.twoD(:);
        else
          zeta2D = zetaStructure.raw.zeta.twoD(:);
        end
            
        median = nanmedian(zeta2D)*1000;
        mean   = nanmean(zeta2D)*1000;
        stdErr = nanserr(zeta2D)*1000;
        stdDev = nanstd(zeta2D)*1000;
        n      = min(sum(~isnan(zeta2D)));

        fprintf('\n2D zeta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'zetaXYZ'
        
        if opts.filtered
          zetaXYZ(:,1) = zetaStructure.filtered.zeta.x(:);
          zetaXYZ(:,2) = zetaStructure.filtered.zeta.y(:);
          zetaXYZ(:,3) = zetaStructure.filtered.zeta.z(:);
        else
          zetaXYZ(:,1) = zetaStructure.raw.zeta.x(:);
          zetaXYZ(:,2) = zetaStructure.raw.zeta.y(:);
          zetaXYZ(:,3) = zetaStructure.raw.zeta.z(:);
        end
        
        median = nanmedian(zetaXYZ)*1000;
        mean   = nanmean(zetaXYZ)*1000;
        stdErr = nanserr(zetaXYZ)*1000;
        stdDev = nanstd(zetaXYZ)*1000;
        n      = min(sum(~isnan(zetaXYZ)));

        fprintf('\nX, Y and Z zeta measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f, %.2f] nm\n',median(1),median(2),median(3))
        fprintf('    mean    = [%.2f, %.2f, %.2f] nm\n',mean(1),mean(2),mean(3))
        fprintf('    std err = [%.2f, %.2f, %.2f] nm\n',stdErr(1),stdErr(2),stdErr(3))
        fprintf('    std dev = [%.2f, %.2f, %.2f] nm\n',stdDev(1),stdDev(2),stdDev(3))
        
    otherwise
        
        error('Statistic requested is not yet built into this version of dublBasicStats. See later release.')
    
end

fprintf('\n')

end