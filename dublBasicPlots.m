function dublBasicPlots(intraStructures,varargin)
% DUBLBASICPLOTS Produces plots of a given statistic of intra-kinetochore
% measurements.
%
%    DUBLBASICPLOTS(INTRASTRUCTURE,...) Outputs a figure representing the
%    distribution of various measurements of intra-kinetochore distances
%    and/or angles, collated in each cell element of INTRASTRUCTURES,
%    across multiple experiments. INTRASTRUCTURES are produced by
%    dublIntraMeasurements.
%
%    Options, defaults in {}:-
%
%    coordSystem: {'plate'} or 'microscope'. Which coordinate system within
%       which to represent measurements. Note that some measurements are
%       only available using 'plate' (e.g. twist), and some measurements
%       are independent of coordinate system (e.g. 3D delta).
%
%    depthFilter: 0 or {1}. Whether or not to give depth-filtered
%       measurements of intra-measurements.
%
%    legend: {'Expt 1', ...} or similar. Names for each experiment.
%    
%    stat: {'delta3D'} or one of the following:
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
% Copyright (c) 2017 C. A. Smith


% default options
opts.coordSystem = 'microscope';
opts.depthFilter = 1;
opts.legend = {'Expt 1','Expt 2','Expt 3','Expt 4', 'Expt 5'};
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

% process input
if ~iscell(intraStructures)
  nExpts = 1;
  intraStructures = {intraStructures};
else
  nExpts = length(intraStructures);
end

% ask the user which stat they would like
statsList = {'delta3D' ,'delta2D' ,'delta1D'  ,'deltaXYZ',...
             'delta3D-church','delta2D-church',...
             'sisSep3D','sisSep2D','sisSepXYZ',...
             'twist3D' ,'twistYZ' ,...
             'swivel3D','swivelYZ','swivelKMT',...
             'rawInts','normInts'};
if ~strcmp(opts.coordSystem,'plate')
    statsList(10:11) = [];
end
if ~ismember(opts.stat,statsList)
    fprintf('Output statistics options:\n')
    for iStat=1:length(statsList)
        fprintf('    %i) %s\n',iStat,statsList{iStat});
    end
    prompt = sprintf('Please type the number statistic you would like: ');
    result = input(prompt);
    opts.stat = statsList{result};
end

for iExpt = 1:nExpts
    isall{iExpt} = isfield(intraStructures{iExpt}.microscope.raw.delta.threeD,'all');
    if strcmp(opts.coordSystem,'plate')
        if (isall{iExpt} && isempty(intraStructures{iExpt}.plate.raw.delta.threeD.all)) || (~isall{iExpt} && isempty(intraStructures{iExpt}.plate.raw.delta.threeD))
            kitLog('No plane fit in movies in intraMeasurements. Converting coordinate system to ''microscope''');
            opts.coordSystem = 'microscope';
        end
    end
end

switch opts.stat

    case 'delta3D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta3D{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.threeD.all(:)*1000;
                        else
                            delta3D{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.threeD(:)*1000;
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta3D{iExpt} = intraStructures{iExpt}.plate.raw.delta.threeD.all(:)*1000;
                        else
                            delta3D{iExpt} = intraStructures{iExpt}.plate.raw.delta.threeD(:)*1000;
                        end
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta3D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.threeD.all(:)*1000;
                        else
                            delta3D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.threeD(:)*1000;
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta3D{iExpt} = intraStructures{iExpt}.microscope.raw.delta.threeD.all(:)*1000;
                        else
                            delta3D{iExpt} = intraStructures{iExpt}.microscope.raw.delta.threeD(:)*1000;
                        end
                    end
                end
        end
        
        title = sprintf('\\Delta_{3D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(delta3D,'axisLimits',[0 300],'title',title,'xLabel','\Delta_{3D} (nm)');
        
    case 'delta2D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta2D{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.twoD.all(:)*1000;
                        else
                            delta2D{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.twoD(:)*1000;
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta2D{iExpt} = intraStructures{iExpt}.plate.raw.delta.twoD.all(:)*1000;
                        else
                            delta2D{iExpt} = intraStructures{iExpt}.plate.raw.delta.twoD(:)*1000;
                        end
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta2D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.twoD.all(:)*1000;
                        else
                            delta2D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.twoD(:)*1000;
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            delta2D{iExpt} = intraStructures{iExpt}.microscope.raw.delta.twoD.all(:)*1000;
                        else
                            delta2D{iExpt} = intraStructures{iExpt}.microscope.raw.delta.twoD(:)*1000;
                        end
                    end
                end
        end

        title = sprintf('\\Delta_{2D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(delta2D,'axisLimits',[0 300],'title',title,'xLabel','\Delta_{2D} (nm)');
        
    case 'delta1D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                      delta1D{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.oneD(:)*1000;
                    end
                else
                    for iExpt = 1:nExpts
                      delta1D{iExpt} = intraStructures{iExpt}.plate.raw.delta.oneD(:)*1000;
                    end
                end   
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                      delta1D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.oneD(:)*1000;
                    end
                else
                    for iExpt = 1:nExpts
                      delta1D{iExpt} = intraStructures{iExpt}.microscope.raw.delta.oneD(:)*1000;
                    end
                end
        end

        title = sprintf('\\Delta_{1D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(delta1D,'axisLimits',[-100 200],'title',title,'xLabel','\Delta_{1D} (nm)');
        
    case 'deltaXYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            deltaX{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.x.all(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.y.all(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.z.all(:)*1000;
                        else
                            deltaX{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.x(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.y(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.plate.depthFilter.delta.z(:)*1000;
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            deltaX{iExpt} = intraStructures{iExpt}.plate.raw.delta.x.all(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.plate.raw.delta.y.all(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.plate.raw.delta.z.all(:)*1000;
                        else
                            deltaX{iExpt} = intraStructures{iExpt}.plate.raw.delta.x(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.plate.raw.delta.y(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.plate.raw.delta.z(:)*1000;
                        end
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            deltaX{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.x.all(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.y.all(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.z.all(:)*1000;
                        else
                            deltaX{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.x(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.y(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.microscope.depthFilter.delta.z(:)*1000;
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            deltaX{iExpt} = intraStructures{iExpt}.microscope.raw.delta.x.all(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.microscope.raw.delta.y.all(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.microscope.raw.delta.z.all(:)*1000;
                        else
                            deltaX{iExpt} = intraStructures{iExpt}.microscope.raw.delta.x(:)*1000;
                            deltaY{iExpt} = intraStructures{iExpt}.microscope.raw.delta.y(:)*1000;
                            deltaZ{iExpt} = intraStructures{iExpt}.microscope.raw.delta.z(:)*1000;
                        end
                    end
                end
        end
        
        figure; clf
        title = sprintf('\\Delta_{x}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        subplot(1,3,1)
        compareHistograms(deltaX,'nBins',10,'axisLimits',[-200 200],...
            'title',title,'xLabel','\Delta_{x} (nm)','withinFig',1);
        subplot(1,3,2)
        title(9) = 'y';
        compareHistograms(deltaY,'nBins',10,'axisLimits',[-200 200],...
            'title',title,'xLabel','\Delta_{y} (nm)','yLabel','','withinFig',1);
        subplot(1,3,3)
        title(9) = 'z';
        compareHistograms(deltaZ,'nBins',10,'axisLimits',[-200 200],...
            'title',title,'xLabel','\Delta_{z} (nm)','yLabel','','withinFig',1);
    
    case 'delta3D-church'
        
        if nExpts > 1
            error('Plotting Churchman-derived inflation-corrected delta only currently works for single conditions.');
        else
            isall = isall{1};
            intraStructure = intraStructures{1};
        end
        
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
        MLp3D(delta3D*1000,[mean stdDev],1,25);
        
    case 'delta2D-church'
        
        if nExpts > 1
            error('Plotting Churchman-derived inflation-corrected delta only currently works for single conditions.');
        else
            isall = isall{1};
            intraStructure = intraStructures{1};
        end
        
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
        MLp2D(delta2D*1000,[mean stdDev],1,25);
        
        
    case 'sisSep3D'
        
        switch opts.coordSystem
            case 'plate'
                for iExpt = 1:nExpts
                  sisSep3D{iExpt} = intraStructures{iExpt}.plate.sisSep.threeD(:);
                end
            case 'microscope'
                for iExpt = 1:nExpts
                  sisSep3D{iExpt} = intraStructures{iExpt}.microscope.sisSep.threeD(:);
                end
        end
        
        title = sprintf('d_{3D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(sisSep3D,'axisLimits',[0.25 2],'title',title,'xLabel','d_{3D} (\mum)');
        
    case 'sisSep2D'
            
        switch opts.coordSystem
            case 'plate'
                for iExpt = 1:nExpts
                  sisSep2D{iExpt} = intraStructures{iExpt}.plate.sisSep.twoD(:);
                end
            case 'microscope'
                for iExpt = 1:nExpts
                  sisSep2D{iExpt} = intraStructures{iExpt}.microscope.sisSep.twoD(:);
                end
        end
        
        title = sprintf('d_{2D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(sisSep2D,'axisLimits',[0.25 1.75],'title',title,'xLabel','d_{2D} (\mum)');
        
    case 'sisSepXYZ'
        
        switch opts.coordSystem
            case 'plate'
                for iExpt = 1:nExpts
                  sisSepX{iExpt} = intraStructures{iExpt}.plate.sisSep.x(:);
                  sisSepY{iExpt} = intraStructures{iExpt}.plate.sisSep.y(:);
                  sisSepZ{iExpt} = intraStructures{iExpt}.plate.sisSep.z(:);
                  xAxisLims = [-3 3];
                end
            case 'microscope'
                for iExpt = 1:nExpts
                  sisSepX{iExpt} = intraStructures{iExpt}.microscope.sisSep.x(:);
                  sisSepY{iExpt} = intraStructures{iExpt}.microscope.sisSep.y(:);
                  sisSepZ{iExpt} = intraStructures{iExpt}.microscope.sisSep.z(:);
                  xAxisLims = [-1.5 1.5];
                end
        end

        figure; clf
        title = sprintf('d_{x}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        subplot(1,3,1)
        compareHistograms(sisSepX,'nBins',10,'axisLimits',xAxisLims,...
            'title',title,'xLabel','d_{x} (nm)','withinFig',1);
        subplot(1,3,2)
        title(4) = 'y';
        compareHistograms(sisSepY,'nBins',10,'axisLimits',[-1.5 1.5],...
            'title',title,'xLabel','d_{y} (nm)','yLabel','','withinFig',1);
        subplot(1,3,3)
        title(4) = 'z';
        compareHistograms(sisSepZ,'nBins',10,'axisLimits',[-1.5 1.5],...
            'title',title,'xLabel','d_{z} (nm)','yLabel','','withinFig',1);

    case 'twist3D'
        
        for iExpt = 1:nExpts
          twist3D{iExpt} = intraStructures{iExpt}.plate.raw.twist.threeD(:);
        end

        figure; clf
        title = sprintf('twist_{3D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(twist3D,'nBins',10,'axisLimits',[0 90],...
            'title',title,'xLabel','twist_{3D} (deg)','yLabel','','withinFig',1);
        
    case 'twistYZ'
        
        for iExpt = 1:nExpts
          twistY{iExpt} = intraStructures{iExpt}.plate.raw.twist.y(:);
          twistZ{iExpt} = intraStructures{iExpt}.plate.raw.twist.z(:);
        end
        
        figure; clf
        title = sprintf('twist_{y}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        subplot(1,2,1)
        compareRoses(twistY,'nBins',18,...
            'title',title,'xLabel','twist_{y} (deg)','withinFig',1);
        subplot(1,2,2)
        title(8) = 'z';
        compareRoses(twistZ,'nBins',18,...
            'title',title,'xLabel','twist_{z} (deg)','withinFig',1);
        
    case 'swivel3D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivel3D{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.threeD.all(:);
                        else
                            swivel3D{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.threeD(:);
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivel3D{iExpt} = intraStructures{iExpt}.plate.raw.swivel.threeD.all(:);
                        else
                            swivel3D{iExpt} = intraStructures{iExpt}.plate.raw.swivel.threeD(:);
                        end
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivel3D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.threeD.all(:);
                        else
                            swivel3D{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.threeD(:);
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivel3D{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.threeD.all(:);
                        else
                            swivel3D{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.threeD(:);
                        end
                    end
                end
        end

        figure; clf
        title = sprintf('swivel_{3D}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareHistograms(swivel3D,'nBins',10,'axisLimits',[0 180],...
            'title',title,'xLabel','swivel_{3D} (deg)','yLabel','','withinFig',1);
        
    case 'swivelYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelY{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.y.all(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.z.all(:);
                        else
                            swivelY{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.y(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.z(:);
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelY{iExpt} = intraStructures{iExpt}.plate.raw.swivel.y.all(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.plate.raw.swivel.z.all(:);
                        else
                            swivelY{iExpt} = intraStructures{iExpt}.plate.raw.swivel.y(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.plate.raw.swivel.z(:);
                        end
                    end
                end        
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelY{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.y.all(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.z.all(:);
                        else
                            swivelY{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.y(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.z(:);
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelY{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.y.all(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.z.all(:);
                        else
                            swivelY{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.y(:);
                            swivelZ{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.z(:);
                        end
                    end
                end  
        end
        figure; clf
        title = sprintf('swivel_{y}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        subplot(1,2,1)
        compareRoses(swivelY,'nBins',18,...
            'title',title,'xLabel','swivel_{y} (deg)','withinFig',1);
        subplot(1,2,2)
        title(9) = 'z';
        compareRoses(swivelZ,'nBins',18,...
            'title',title,'xLabel','swivel_{z} (deg)','withinFig',1);
        
    case 'swivelKMT'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelKMT{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.kMT.all(:);
                        else
                            swivelKMT{iExpt} = intraStructures{iExpt}.plate.depthFilter.swivel.kMT(:);
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelKMT{iExpt} = intraStructures{iExpt}.plate.raw.swivel.kMT.all(:);
                        else
                            swivelKMT{iExpt} = intraStructures{iExpt}.plate.raw.swivel.kMT(:);
                        end
                    end
                end        
            case 'microscope'
                if opts.depthFilter
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelKMT{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.kMT.all(:);
                        else
                            swivelKMT{iExpt} = intraStructures{iExpt}.microscope.depthFilter.swivel.kMT(:);
                        end
                    end
                else
                    for iExpt = 1:nExpts
                        if isall{iExpt}
                            swivelKMT{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.kMT.all(:);
                        else
                            swivelKMT{iExpt} = intraStructures{iExpt}.microscope.raw.swivel.kMT(:);
                        end
                    end
                end  
        end
        figure; clf
        title = sprintf('swivel_{kMT}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareRoses(swivelY,'nBins',18,...
            'title',title,'xLabel','swivel_{kMT} (deg)','withinFig',1);
        
    case 'rawInts'
        
        for iExpt = 1:nExpts
            intsInner{iExpt} = intraStructures{iExpt}.intensity.mean.inner(:);
            intsOuter{iExpt} = intraStructures{iExpt}.intensity.mean.outer(:);
            if opts.depthFilter
                if isall{iExpt}
                    filt = ~isnan(intraStructures{iExpt}.plate.depthFilter.delta.threeD.all(:));
                else
                    filt = ~isnan(intraStructures{iExpt}.plate.depthFilter.delta.threeD(:));
                end
                intsInner{iExpt}(filt) = NaN;
                intsOuter{iExpt}(filt) = NaN;
            end
        end
        
        figure; clf
        subplot(1,2,1);
        title = sprintf('raw intensity_{inner}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareBoxWhiskers(intsInner,'showOutliers',0,'outlierP',0.05,...
            'title',title,'yLabel','raw intensity',...
            'withinFig',1,'legend',opts.legend);
        subplot(1,2,2);
        title = sprintf('raw intensity_{outer}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareBoxWhiskers(intsOuter,'showOutliers',0,'outlierP',0.05,...
            'title',title,'yLabel','raw intensity',...
            'withinFig',1,'legend',opts.legend);
        
    case 'normInts'
        
        for iExpt = 1:nExpts
            intsInner{iExpt} = intraStructures{iExpt}.intensity.mean.inner(:);
            intsOuter{iExpt} = intraStructures{iExpt}.intensity.mean.outer(:);
            if opts.depthFilter
                if isall{iExpt}
                    filt = ~isnan(intraStructures{iExpt}.plate.depthFilter.delta.threeD.all(:));
                else
                    filt = ~isnan(intraStructures{iExpt}.plate.depthFilter.delta.threeD(:));
                end
                intsInner{iExpt}(filt) = NaN;
                intsOuter{iExpt}(filt) = NaN;
            end
            normIntsInner{iExpt} = intsInner{iExpt};
            normIntsOuter{iExpt} = intsOuter{iExpt};
            normIntsInner{iExpt} = normIntsInner{iExpt}./intsOuter{iExpt};
            normIntsOuter{iExpt} = normIntsOuter{iExpt}./intsInner{iExpt};
        end
        
        figure; clf
        subplot(1,2,1);
        title = sprintf('inner-normalised intensity_{outer}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareBoxWhiskers(normIntsOuter,'showOutliers',0,'outlierP',0.05,...
            'title',title,'yLabel','normalised intensity',...
            'withinFig',1,'legend',opts.legend);
        subplot(1,2,2);
        title = sprintf('outer-normalised intensity_{inner}: %s',opts.legend{1});
        for iExpt = 2:nExpts
          title = [title ' vs. ' opts.legend{iExpt}];
        end
        compareBoxWhiskers(normIntsInner,'showOutliers',0,'outlierP',0.05,...
            'title',title,'yLabel','normalised intensity',...
            'withinFig',1,'legend',opts.legend);
        
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