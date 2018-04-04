function dublBasicStats(intraStructure,varargin)
% DUBLBASICSTATS Gives a list of statistics for delta measurements.
%
%    DUBLBASICSTATS(COMPILEDDELTA,...) Calculates the median, mean and
%    standard error and deviation of measurements of delta collated in
%    INTRASTRUCTURE and outputs to the command line. INTRASTRUCTURE is
%    produced by dublIntraMeasurements.
%
%    Options, defaults in {}:-
%
%    coordSystem: {'plate'} or 'microscope'. The coordinate system in which
%       to provide statistics. Note that some statistics are coordinate-
%       independent.
%
%    depthFilter: 0 or {1}. Whether or not to give depth-filtered
%       measurements of intra-measurements.
%
%    stat: {''} or one of the following:
%           - 'delta3D'
%           - 'delta2D'
%           - 'delta1D'
%           - 'deltaXYZ'
%           - 'delta3D-church'
%           - 'delta2D-church'
%           - 'sisSep3D'
%           - 'sisSep2D'
%           - 'sisSepXYZ'
%           - 'twist3D'
%           - 'twistYZ'
%           - 'swivel3D'
%           - 'swivelYZ'
%           - 'swivelKMT'
%           - 'rawInts'
%           - 'normInts'
%       The statistic to be printed to screen. If no statistic is provided,
%       the user will be prompted.
%
% Copyright (c) 2018 C. A. Smith


% default options
opts.coordSystem = 'microscope';
opts.depthFilter = 1;
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

statsList = {'delta3D' ,'delta2D' ,'delta1D'  ,'deltaXYZ',...
             'delta3D-church','delta2D-church',...
             'sisSep3D','sisSep2D','sisSepXYZ',...
             'twist3D' ,'twistYZ' ,...
             'swivel3D','swivelYZ','swivelKMT',...
             'rawInts','normInts'};
if ~ismember(opts.stat,statsList)
    fprintf('Output statistics options:\n')
    for iStat=1:length(statsList)
        fprintf('    %i) %s\n',iStat,statsList{iStat});
    end
    prompt = sprintf('Please type the number statistic you would like: ');
    result = input(prompt);
    opts.stat = statsList{result};
end
isall = isfield(intraStructure.microscope.raw.delta.threeD,'all');

if strcmp(opts.coordSystem,'plate')
    if (isall && isempty(intraStructure.plate.raw.delta.threeD.all)) || (~isall && isempty(intraStructure.plate.raw.delta.threeD))
        kitLog('No plane fit in movies in intraMeasurements. Converting coordinate system to ''microscope''');
        opts.coordSystem = 'microscope';
    end
end

switch opts.stat

    case 'delta3D'
            
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        delta3D = intraStructure.plate.depthFilter.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.plate.depthFilter.delta.threeD(:);
                    end
                else
                    if isall
                        delta3D = intraStructure.plate.raw.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.plate.raw.delta.threeD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        delta3D = intraStructure.microscope.depthFilter.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.microscope.depthFilter.delta.threeD(:);
                    end
                else
                    if isall
                        delta3D = intraStructure.microscope.raw.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.microscope.raw.delta.threeD(:);
                    end
                end
        end
        
        median = nanmedian(delta3D)*1000;
        mean   = nanmean(delta3D)*1000;
        stdErr = nanserr(delta3D)*1000;
        stdDev = nanstd(delta3D)*1000;
        n      = min(sum(~isnan(delta3D)));

        fprintf('\n3D delta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'delta2D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        delta2D = intraStructure.plate.depthFilter.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.plate.depthFilter.delta.twoD(:);
                    end
                else
                    if isall
                        delta2D = intraStructure.plate.raw.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.plate.raw.delta.twoD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        delta2D = intraStructure.microscope.depthFilter.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.microscope.depthFilter.delta.twoD(:);
                    end
                else
                    if isall
                        delta2D = intraStructure.microscope.raw.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.microscope.raw.delta.twoD(:);
                    end
                end
        end

        median = nanmedian(delta2D)*1000;
        mean   = nanmean(delta2D)*1000;
        stdErr = nanserr(delta2D)*1000;
        stdDev = nanstd(delta2D)*1000;
        n      = min(sum(~isnan(delta2D)));

        fprintf('\n2D delta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'delta1D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    delta1D = intraStructure.plate.depthFilter.delta.oneD(:);
                else
                    delta1D = intraStructure.plate.raw.delta.oneD(:);
                end   
            case 'microscope'
                if opts.depthFilter
                    delta1D = intraStructure.microscope.depthFilter.delta.oneD(:);
                else
                    delta1D = intraStructure.microscope.raw.delta.oneD(:);
                end
        end

        median = nanmedian(delta1D)*1000;
        mean   = nanmean(delta1D)*1000;
        stdErr = nanserr(delta1D)*1000;
        stdDev = nanstd(delta1D)*1000;
        n      = min(sum(~isnan(delta1D)));

        fprintf('\n1D delta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'deltaXYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z(:);
                    end
                else
                    if isall
                        deltaXYZ(:,1) = intraStructure.plate.raw.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.plate.raw.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.plate.raw.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.plate.raw.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.plate.raw.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.plate.raw.delta.z(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z(:);
                    end
                else
                    if isall
                        deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z(:);
                    end
                end
        end
    
        median = nanmedian(deltaXYZ)*1000;
        mean   = nanmean(deltaXYZ)*1000;
        stdErr = nanserr(deltaXYZ)*1000;
        stdDev = nanstd(deltaXYZ)*1000;
        n      = min(sum(~isnan(deltaXYZ)));

        fprintf('\nX, Y and Z delta measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f, %.2f] nm\n',median(1),median(2),median(3))
        fprintf('    mean    = [%.2f, %.2f, %.2f] nm\n',mean(1),mean(2),mean(3))
        fprintf('    std err = [%.2f, %.2f, %.2f] nm\n',stdErr(1),stdErr(2),stdErr(3))
        fprintf('    std dev = [%.2f, %.2f, %.2f] nm\n',stdDev(1),stdDev(2),stdDev(3))
    
    case 'delta3D-church'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        delta3D = intraStructure.plate.depthFilter.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.plate.depthFilter.delta.threeD(:);
                    end
                else
                    if isall
                        delta3D = intraStructure.plate.raw.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.plate.raw.delta.threeD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        delta3D = intraStructure.microscope.depthFilter.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.microscope.depthFilter.delta.threeD(:);
                    end
                else
                    if isall
                        delta3D = intraStructure.microscope.raw.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.microscope.raw.delta.threeD(:);
                    end
                end
        end
        
        % Churchman falls over with outliers - remove them.
        outs = findoutliers(delta3D);
        delta3D = delta3D(~outs);
        
        mean   = nanmean(delta3D)*1000;
        stdDev = nanstd(delta3D)*1000;
        n      = sum(~isnan(delta3D));
        [church,~] = MLp3D(delta3D*1000,[mean stdDev]);
        
        fprintf('\nChurchman-corrected 3D delta measurements (n = %i):\n\n',n)
        fprintf('    mean    = %.2f nm\n',church(1))
        fprintf('    std dev = %.2f nm\n',church(2)) 
        
    case 'delta2D-church'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        delta2D = intraStructure.plate.depthFilter.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.plate.depthFilter.delta.twoD(:);
                    end
                else
                    if isall
                        delta2D = intraStructure.plate.raw.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.plate.raw.delta.twoD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        delta2D = intraStructure.microscope.depthFilter.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.microscope.depthFilter.delta.twoD(:);
                    end
                else
                    if isall
                        delta2D = intraStructure.microscope.raw.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.microscope.raw.delta.twoD(:);
                    end
                end
        end
        
        % Churchman falls over with outliers - remove them.
        outs = findoutliers(delta2D);
        delta2D = delta2D(~outs);
        
        mean   = nanmean(delta2D)*1000;
        stdDev = nanstd(delta2D)*1000;
        n      = sum(~isnan(delta2D));
        [church,~] = MLp2D(delta2D*1000,[mean stdDev]);
        
        fprintf('\nChurchman-corrected 2D delta measurements (n = %i):\n\n',n)
        fprintf('    mean    = %.2f nm\n',church(1))
        fprintf('    std dev = %.2f nm\n',church(2)) 
        
    case 'sisSep3D'
        
        switch opts.coordSystem
            case 'plate'
                sisSep3D = intraStructure.plate.sisSep.threeD(:);
            case 'microscope'
                sisSep3D = intraStructure.microscope.sisSep.threeD(:);
        end
        
        median = nanmedian(sisSep3D);
        mean   = nanmean(sisSep3D);
        stdErr = nanserr(sisSep3D);
        stdDev = nanstd(sisSep3D);
        n      = min(sum(~isnan(sisSep3D)));

        fprintf('\n3D sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f um\n',median)
        fprintf('    mean    = %.2f um\n',mean)
        fprintf('    std err = %.2f um\n',stdErr)
        fprintf('    std dev = %.2f um\n',stdDev)
        
    case 'sisSep2D'
            
        switch opts.coordSystem
            case 'plate'
                sisSep2D = intraStructure.plate.sisSep.twoD(:);
            case 'microscope'
                sisSep2D = intraStructure.microscope.sisSep.twoD(:);
        end
        
        median = nanmedian(sisSep2D);
        mean   = nanmean(sisSep2D);
        stdErr = nanserr(sisSep2D);
        stdDev = nanstd(sisSep2D);
        n      = min(sum(~isnan(sisSep2D)));

        fprintf('\n2D sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f um\n',median)
        fprintf('    mean    = %.2f um\n',mean)
        fprintf('    std err = %.2f um\n',stdErr)
        fprintf('    std dev = %.2f um\n',stdDev)
        
    case 'sisSepXYZ'
        
        switch opts.coordSystem
            case 'plate'
                sisSepXYZ(:,1) = intraStructure.plate.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.plate.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.plate.sisSep.z(:);
            case 'microscope'
                sisSepXYZ(:,1) = intraStructure.microscope.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.microscope.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.microscope.sisSep.z(:);
        end

        median = nanmedian(sisSepXYZ);
        mean   = nanmean(sisSepXYZ);
        stdErr = nanserr(sisSepXYZ);
        stdDev = nanstd(sisSepXYZ);
        n      = min(sum(~isnan(sisSepXYZ)));

        fprintf('\nX, Y and Z sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f, %.2f] um\n',median(1),median(2),median(3))
        fprintf('    mean    = [%.2f, %.2f, %.2f] um\n',mean(1),mean(2),mean(3))
        fprintf('    std err = [%.2f, %.2f, %.2f] um\n',stdErr(1),stdErr(2),stdErr(3))
        fprintf('    std dev = [%.2f, %.2f, %.2f] um\n',stdDev(1),stdDev(2),stdDev(3))

    case 'twist3D'
        
        switch opts.coordSystem
            case 'plate'
                twist3D = intraStructure.plate.twist.threeD(:);
            case 'microscope'
                error('Twist is not defined for a non-plate coordinate system.')
        end
        
        median = nanmedian(twist3D);
        mean   = nanmean(twist3D);
        stdErr = nanserr(twist3D);
        stdDev = nanstd(twist3D);
        n      = min(sum(~isnan(twist3D)));

        fprintf('\n3D twist measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f deg\n',median)
        fprintf('    mean    = %.2f deg\n',mean)
        fprintf('    std err = %.2f deg\n',stdErr)
        fprintf('    std dev = %.2f deg\n',stdDev)
        
    case 'twistYZ'
        
        switch opts.coordSystem
            case 'plate'
                twistYZ(:,1) = intraStructure.plate.twist.y(:);
                twistYZ(:,2) = intraStructure.plate.twist.z(:);
            case 'microscope'
                error('Twist is not defined for a non-plate coordinate system.')
        end

        median = nanmedian(twistYZ);
        mean   = nanmean(twistYZ);
        stdErr = nanserr(twistYZ);
        stdDev = nanstd(twistYZ);
        n      = min(sum(~isnan(twistYZ)));

        fprintf('\nY and Z twist measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] deg\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] deg\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] deg\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] deg\n',stdDev(1),stdDev(2))

    case 'swivel3D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        swivel3D = intraStructure.plate.depthFilter.swivel.threeD.all(:);
                    else
                        swivel3D = intraStructure.plate.depthFilter.swivel.threeD(:);
                    end
                else
                    if isall
                        swivel3D = intraStructure.plate.raw.swivel.threeD.all(:);
                    else
                        swivel3D = intraStructure.plate.raw.swivel.threeD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        swivel3D = intraStructure.microscope.depthFilter.swivel.threeD.all(:);
                    else
                        swivel3D = intraStructure.microscope.depthFilter.swivel.threeD(:);
                    end
                else
                    if isall
                        swivel3D = intraStructure.microscope.raw.swivel.threeD.all(:);
                    else
                        swivel3D = intraStructure.microscope.raw.swivel.threeD(:);
                    end
                end
        end
        
        median = nanmedian(swivel3D);
        mean   = nanmean(swivel3D);
        stdErr = nanserr(swivel3D);
        stdDev = nanstd(swivel3D);
        n      = min(sum(~isnan(swivel3D)));

        fprintf('\n3D swivel measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f deg\n',median)
        fprintf('    mean    = %.2f deg\n',mean)
        fprintf('    std err = %.2f deg\n',stdErr)
        fprintf('    std dev = %.2f deg\n',stdDev)
        
    case 'swivelKMT'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    swivelKMT = intraStructure.plate.depthFilter.swivel.kMT(:);
                else
                    swivelKMT = intraStructure.plate.raw.swivel.kMT(:);
                end
            case 'microscope'
                if opts.depthFilter
                    swivelKMT = intraStructure.microscope.depthFilter.swivel.kMT(:);
                else
                    swivelKMT = intraStructure.microscope.raw.swivel.kMT(:);
                end
        end

        median = nanmedian(swivelKMT);
        mean   = nanmean(swivelKMT);
        stdErr = nanserr(swivelKMT);
        stdDev = nanstd(swivelKMT);
        n      = min(sum(~isnan(swivelKMT)));

        fprintf('\nKinetochore-to-microtubule angle measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f deg\n',median)
        fprintf('    mean    = %.2f deg\n',mean)
        fprintf('    std err = %.2f deg\n',stdErr)
        fprintf('    std dev = %.2f deg\n',stdDev)
        
    case 'swivelYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        swivelYZ(:,1) = intraStructure.plate.depthFilter.swivel.y.all(:);
                        swivelYZ(:,2) = intraStructure.plate.depthFilter.swivel.z.all(:);
                    else
                        swivelYZ(:,1) = intraStructure.plate.depthFilter.swivel.y(:);
                        swivelYZ(:,2) = intraStructure.plate.depthFilter.swivel.z(:);
                    end
                else
                    if isall
                        swivelYZ(:,1) = intraStructure.plate.raw.swivel.y.all(:);
                        swivelYZ(:,2) = intraStructure.plate.raw.swivel.z.all(:);
                    else
                        swivelYZ(:,1) = intraStructure.plate.raw.swivel.y(:);
                        swivelYZ(:,2) = intraStructure.plate.raw.swivel.z(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        swivelYZ(:,1) = intraStructure.microscope.depthFilter.swivel.y.all(:);
                        swivelYZ(:,2) = intraStructure.microscope.depthFilter.swivel.z.all(:);
                    else
                        swivelYZ(:,1) = intraStructure.microscope.depthFilter.swivel.y(:);
                        swivelYZ(:,2) = intraStructure.microscope.depthFilter.swivel.z(:);
                    end
                else
                    if isall
                        swivelYZ(:,1) = intraStructure.microscope.raw.swivel.y.all(:);
                        swivelYZ(:,2) = intraStructure.microscope.raw.swivel.z.all(:);
                    else
                        swivelYZ(:,1) = intraStructure.microscope.raw.swivel.y(:);
                        swivelYZ(:,2) = intraStructure.microscope.raw.swivel.z(:);
                    end
                end
        end
    
        median = nanmedian(swivelYZ);
        mean   = nanmean(swivelYZ);
        stdErr = nanserr(swivelYZ);
        stdDev = nanstd(swivelYZ);
        n      = min(sum(~isnan(swivelYZ)));

        fprintf('\nY and Z swivel measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] deg\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] deg\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] deg\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] deg\n',stdDev(1),stdDev(2))
        
    case 'rawInts'
        
        ints(:,1) = intraStructure.intensity.mean.inner(:);
        ints(:,2) = intraStructure.intensity.mean.outer(:);
        if opts.depthFilter
            if isall
                filt = ~isnan(intraStructure.plate.depthFilter.delta.threeD.all(:));
            else
                filt = ~isnan(intraStructure.plate.depthFilter.delta.threeD(:));
            end
            ints(repmat(filt,1,2)) = NaN;
        end
        
        median = nanmedian(ints)*1000;
        mean   = nanmean(ints)*1000;
        stdErr = nanserr(ints)*1000;
        stdDev = nanstd(ints)*1000;
        n      = min(sum(~isnan(ints)));

        fprintf('\nRaw inner- and outer-marker intensity measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.4f, %.4f] au\n',median(1),median(2))
        fprintf('    mean    = [%.4f, %.4f] au\n',mean(1),mean(2))
        fprintf('    std err = [%.4f, %.4f] au\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.4f, %.4f] au\n',stdDev(1),stdDev(2))
        
    case 'normInts'
        
        ints(:,1) = intraStructure.intensity.mean.inner(:);
        ints(:,2) = intraStructure.intensity.mean.outer(:);
        if opts.depthFilter
            if isall
                filt = ~isnan(intraStructure.plate.depthFilter.delta.threeD.all(:));
            else
                filt = ~isnan(intraStructure.plate.depthFilter.delta.threeD(:));
            end
            ints(repmat(filt,1,2)) = NaN;
        end
        normInts = ints;
        normInts(:,2) = normInts(:,2)./ints(:,1);
        normInts(:,1) = normInts(:,1)./ints(:,2);
        
        median = nanmedian(normInts);
        mean   = nanmean(normInts);
        stdErr = nanserr(normInts);
        stdDev = nanstd(normInts);
        n      = min(sum(~isnan(normInts)));

        fprintf('\nInner-normalised outer-marker intensity measurements (n = %i):\n\n',n)
        fprintf('    median  = %.4f au\n',median(2))
        fprintf('    mean    = %.4f au\n',mean(2))
        fprintf('    std err = %.4f au\n',stdErr(2))
        fprintf('    std dev = %.4f au\n',stdDev(2))
        fprintf('\n');
        fprintf('\nOuter-normalised inner-marker intensity measurements (n = %i):\n\n',n)
        fprintf('    median  = %.4f au\n',median(1))
        fprintf('    mean    = %.4f au\n',mean(1))
        fprintf('    std err = %.4f au\n',stdErr(1))
        fprintf('    std dev = %.4f au\n',stdDev(1))
        
    otherwise
        
        error('Statistic requested is not yet built into this version of dublBasicStats. See later release.')
    
end

fprintf('\n')

end

%% Sub-functions

function outs = findoutliers(data)
  if nargin<1 || isempty(data)
    return
  end
  if verLessThan('matlab','9.2')
    nTests = length(data);
    outs = zeros(nTests,1);
    for iTest = 1:nTests
      outs(iTest) = ttest2(data,data(iTest),'alpha',0.0455);
    end    
  else
    outs = isoutlier(data,'mean');
  end
end