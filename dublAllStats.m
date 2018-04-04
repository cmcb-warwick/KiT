function dublAllStats(intraStructure,varargin)
% DUBLALLSTATS Gives a list of statistics for delta measurements.
%
%    DUBLALLSTATS(INTRASTRUCTURE,...) Calculates the median, mean, standard
%    error and deviation, and number of measurements collated in
%    INTRASTRUCTURE and outputs a table to the command line. INTRASTRUCTURE
%    is produced by dublIntraMeasurements.
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
%
% Copyright (c) 2018 C. A. Smith


% default options
opts.coordSystem = 'microscope';
opts.depthFilter = 1;
% get user options
opts = processOptions(opts,varargin{:});

ispaired = isfield(intraStructure.microscope.raw,'swivel');
if ispaired
    statsList = {'delta3D';'delta2D';'delta1D';'deltaX';'deltaY';'deltaZ';...
             'delta3D-church';'delta2D-church';...
             'sisSep3D';'sisSep2D';'sisSepX';'sisSepY';'sisSepZ';...
             'twist3D';'twistY';'twistZ';...
             'swivel3D';'swivelY';'swivelZ';'swivelKMT'};
    Median = nan(20,1);
    Mean   = nan(20,1);
    StdErr = nan(20,1);
    StdDev = nan(20,1);
    n      = nan(20,1);
else
    statsList = {'delta3D';'delta2D';'deltaX';'deltaY';'deltaZ';...
             'delta3D-church';'delta2D-church'};
    Median = nan(7,1);
    Mean   = nan(7,1);
    StdErr = nan(7,1);
    StdDev = nan(7,1);
    n      = nan(7,1);
end
         
isall = isfield(intraStructure.microscope.raw.delta.threeD,'all');

if strcmp(opts.coordSystem,'plate')
    if (isall && isempty(intraStructure.plate.raw.delta.threeD.all)) || (~isall && isempty(intraStructure.plate.raw.delta.threeD))
        kitLog('No plane fit in movies in intraMeasurements. Converting coordinate system to ''microscope''');
        opts.coordSystem = 'microscope';
    end
end

switch opts.coordSystem
    case 'plate'
        if opts.depthFilter
            if isall
                delta3D = intraStructure.plate.depthFilter.delta.threeD.all(:);
                delta2D = intraStructure.plate.depthFilter.delta.twoD.all(:);
                deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x.all(:);
                deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y.all(:);
                deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z.all(:);
            else
                delta3D = intraStructure.plate.depthFilter.delta.threeD(:);
                delta2D = intraStructure.plate.depthFilter.delta.twoD(:);
                deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x(:);
                deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y(:);
                deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z(:);
            end
            if ispaired
                if isall
                    swivel3D = intraStructure.plate.depthFilter.swivel.threeD.all(:);
                    swivelYZ(:,1) = intraStructure.plate.depthFilter.swivel.y.all(:);
                    swivelYZ(:,2) = intraStructure.plate.depthFilter.swivel.z.all(:);
                else
                    swivel3D = intraStructure.plate.depthFilter.swivel.threeD(:);
                    swivelYZ(:,1) = intraStructure.plate.depthFilter.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.plate.depthFilter.swivel.z(:);
                end
                delta1D = intraStructure.plate.depthFilter.delta.oneD(:);
                swivelKMT = intraStructure.plate.depthFilter.swivel.kMT(:);
            end
        else
            if isall
                delta3D = intraStructure.plate.raw.delta.threeD.all(:);
                delta2D = intraStructure.plate.raw.delta.twoD.all(:);
                deltaXYZ(:,1) = intraStructure.plate.raw.delta.x.all(:);
                deltaXYZ(:,2) = intraStructure.plate.raw.delta.y.all(:);
                deltaXYZ(:,3) = intraStructure.plate.raw.delta.z.all(:);
            else
                delta3D = intraStructure.plate.raw.delta.threeD(:);
                delta2D = intraStructure.plate.raw.delta.twoD(:);
                deltaXYZ(:,1) = intraStructure.plate.raw.delta.x(:);
                deltaXYZ(:,2) = intraStructure.plate.raw.delta.y(:);
                deltaXYZ(:,3) = intraStructure.plate.raw.delta.z(:);
            end
            if ispaired
                if isall
                    swivel3D = intraStructure.plate.raw.swivel.threeD.all(:);
                    swivelYZ(:,1) = intraStructure.plate.raw.swivel.y.all(:);
                    swivelYZ(:,2) = intraStructure.plate.raw.swivel.z.all(:);
                else
                    swivel3D = intraStructure.plate.raw.swivel.threeD(:);
                    swivelYZ(:,1) = intraStructure.plate.raw.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.plate.raw.swivel.z(:);
                end
                delta1D = intraStructure.plate.raw.delta.oneD(:);
                swivelKMT = intraStructure.plate.raw.swivel.kMT(:);
            end
        end
        if ispaired
            sisSep3D = intraStructure.plate.sisSep.threeD(:);
            sisSep2D = intraStructure.plate.sisSep.twoD(:);
            sisSepXYZ(:,1) = intraStructure.plate.sisSep.x(:);
            sisSepXYZ(:,2) = intraStructure.plate.sisSep.y(:);
            sisSepXYZ(:,3) = intraStructure.plate.sisSep.z(:);
            if isfield(intraStructure.plate,'twist')
                twist3D = intraStructure.plate.twist.threeD(:);
                twistYZ(:,1) = intraStructure.plate.twist.y(:);
                twistYZ(:,2) = intraStructure.plate.twist.z(:);
            else
                twist3D = NaN;
                twistYZ = [NaN NaN];
            end
        end
    case 'microscope'
        if opts.depthFilter
            if isall
                delta3D = intraStructure.microscope.depthFilter.delta.threeD.all(:);
                delta2D = intraStructure.microscope.depthFilter.delta.twoD.all(:);
                deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x.all(:);
                deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y.all(:);
                deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z.all(:);
            else
                delta3D = intraStructure.microscope.depthFilter.delta.threeD(:);
                delta2D = intraStructure.microscope.depthFilter.delta.twoD(:);
                deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x(:);
                deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y(:);
                deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z(:);
            end
            if ispaired
                if isall
                    swivel3D = intraStructure.microscope.depthFilter.swivel.threeD.all(:);
                    swivelYZ(:,1) = intraStructure.microscope.depthFilter.swivel.y.all(:);
                    swivelYZ(:,2) = intraStructure.microscope.depthFilter.swivel.z.all(:);
                else
                    swivel3D = intraStructure.microscope.depthFilter.swivel.threeD(:);
                    swivelYZ(:,1) = intraStructure.microscope.depthFilter.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.microscope.depthFilter.swivel.z(:);
                end
                delta1D = intraStructure.microscope.depthFilter.delta.oneD(:);
                swivelKMT = intraStructure.microscope.depthFilter.swivel.kMT(:);
            end
        else
            if isall
                delta3D = intraStructure.microscope.raw.delta.threeD.all(:);
                delta2D = intraStructure.microscope.raw.delta.twoD.all(:);
                deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x.all(:);
                deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y.all(:);
                deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z.all(:);
            else
                delta3D = intraStructure.microscope.raw.delta.threeD(:);
                delta2D = intraStructure.microscope.raw.delta.twoD(:);
                deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x(:);
                deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y(:);
                deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z(:);
            end
            if ispaired
                if isall
                    swivel3D = intraStructure.microscope.raw.swivel.threeD.all(:);
                    swivelYZ(:,1) = intraStructure.microscope.raw.swivel.y.all(:);
                    swivelYZ(:,2) = intraStructure.microscope.raw.swivel.z.all(:);
                else
                    swivel3D = intraStructure.microscope.raw.swivel.threeD(:);
                    swivelYZ(:,1) = intraStructure.microscope.raw.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.microscope.raw.swivel.z(:);
                end
                delta1D = intraStructure.microscope.raw.delta.oneD(:);
                swivelKMT = intraStructure.microscope.raw.swivel.kMT(:);
            end
        end
        if ispaired
            sisSep3D = intraStructure.microscope.sisSep.threeD(:);
            sisSep2D = intraStructure.microscope.sisSep.twoD(:);
            sisSepXYZ(:,1) = intraStructure.microscope.sisSep.x(:);
            sisSepXYZ(:,2) = intraStructure.microscope.sisSep.y(:);
            sisSepXYZ(:,3) = intraStructure.microscope.sisSep.z(:);
            twist3D = NaN;
            twistYZ = [NaN NaN];
        end
end

% Churchman falls over with outliers - remove them.
chdelta3D = delta3D;
outs = findoutliers(chdelta3D);
chdelta3D = chdelta3D(~outs);

mean   = nanmean(chdelta3D)*1000;
stdDev = nanstd(chdelta3D)*1000;
nch3D      = sum(~isnan(chdelta3D));
[church3D,~] = MLp3D(chdelta3D*1000,[mean stdDev]);

% Churchman falls over with outliers - remove them.
chdelta2D = delta2D;
outs = findoutliers(chdelta2D);
chdelta2D = chdelta2D(~outs);

mean   = nanmean(chdelta2D)*1000;
stdDev = nanstd(chdelta2D)*1000;
nch2D      = sum(~isnan(chdelta2D));
[church2D,~] = MLp2D(chdelta2D*1000,[mean stdDev]);

c=0;

c=c+1;
Median(c) = nanmedian(delta3D)*1000;
Mean(c)   = nanmean(delta3D)*1000;
StdErr(c) = nanserr(delta3D)*1000;
StdDev(c) = nanstd(delta3D)*1000;
n(c,:)      = min(sum(~isnan(delta3D)));
    
c=c+1;
Median(c) = nanmedian(delta2D)*1000;
Mean(c)   = nanmean(delta2D)*1000;
StdErr(c) = nanserr(delta2D)*1000;
StdDev(c) = nanstd(delta2D)*1000;
n(c)      = min(sum(~isnan(delta2D)));

if ispaired
    c=c+1;
    Median(c) = nanmedian(delta1D)*1000;
    Mean(c)   = nanmean(delta1D)*1000;
    StdErr(c) = nanserr(delta1D)*1000;
    StdDev(c) = nanstd(delta1D)*1000;
    n(c)      = min(sum(~isnan(delta1D)));
end

c=c+1;
Median(c:c+2) = nanmedian(deltaXYZ)*1000;
Mean(c:c+2)   = nanmean(deltaXYZ)*1000;
StdErr(c:c+2) = nanserr(deltaXYZ)*1000;
StdDev(c:c+2) = nanstd(deltaXYZ)*1000;
n(c:c+2)      = min(sum(~isnan(deltaXYZ)));
c=c+2;

c=c+1;
Median(c) = NaN;
Mean(c)   = church3D(1);
StdErr(c) = NaN;
StdDev(c) = church3D(2);
n(c)      = nch3D;

c=c+1;
Median(c) = NaN;
Mean(c)   = church2D(1);
StdErr(c) = NaN;
StdDev(c) = church2D(2);
n(c)      = nch2D;

if ispaired
    c=c+1;
    Median(c) = nanmedian(sisSep3D);
    Mean(c)   = nanmean(sisSep3D);
    StdErr(c) = nanserr(sisSep3D);
    StdDev(c) = nanstd(sisSep3D);
    n(c)      = min(sum(~isnan(sisSep3D)));

    c=c+1;
    Median(c) = nanmedian(sisSep2D);
    Mean(c)   = nanmean(sisSep2D);
    StdErr(c) = nanserr(sisSep2D);
    StdDev(c) = nanstd(sisSep2D);
    n(c)      = min(sum(~isnan(sisSep2D)));

    c=c+1;
    Median(c:c+2) = nanmedian(sisSepXYZ);
    Mean(c:c+2)   = nanmean(sisSepXYZ);
    StdErr(c:c+2) = nanserr(sisSepXYZ);
    StdDev(c:c+2) = nanstd(sisSepXYZ);
    n(c:c+2)      = min(sum(~isnan(sisSepXYZ)));
    c=c+2;

    c=c+1;
    Median(c) = nanmedian(twist3D);
    Mean(c)   = nanmean(twist3D);
    StdErr(c) = nanserr(twist3D);
    StdDev(c) = nanstd(twist3D);
    n(c)      = min(sum(~isnan(twist3D)));

    c=c+1;
    Median(c:c+1) = nanmedian(twistYZ);
    Mean(c:c+1)   = nanmean(twistYZ);
    StdErr(c:c+1) = nanserr(twistYZ);
    StdDev(c:c+1) = nanstd(twistYZ);
    n(c:c+1)      = min(sum(~isnan(twistYZ)));
    c=c+1;

    c=c+1;
    Median(c) = nanmedian(swivel3D);
    Mean(c)   = nanmean(swivel3D);
    StdErr(c) = nanserr(swivel3D);
    StdDev(c) = nanstd(swivel3D);
    n(c)      = min(sum(~isnan(swivel3D)));

    c=c+1;
    Median(c:c+1) = nanmedian(swivelYZ);
    Mean(c:c+1)   = nanmean(swivelYZ);
    StdErr(c:c+1) = nanserr(swivelYZ);
    StdDev(c:c+1) = nanstd(swivelYZ);
    n(c:c+1)      = min(sum(~isnan(swivelYZ)));
    c=c+1;

    c=c+1;
    Median(c) = nanmedian(swivelKMT);
    Mean(c)   = nanmean(swivelKMT);
    StdErr(c) = nanserr(swivelKMT);
    StdDev(c) = nanstd(swivelKMT);
    n(c)      = min(sum(~isnan(swivelKMT)));
end

allStats = table(Median,Mean,StdDev,StdErr,n,'RowNames',statsList);
fprintf('\n');
disp(allStats);

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