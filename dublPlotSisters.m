function dublPlotSisters(job, varargin)
%DUBLPLOTSISTERS Produces a set of plots as in kitBasicPlots but for both
%              channels.
%
%    DUBLPLOTSISTERS(JOB,...) Plots all sister pairs for two channels
%    overlaid. By default also plots intersister distances and intermarker
%    distances. Supply options as string/value pairs following JOB.
%
%    Options, defaults in {}:-
%
%    channels: {1:2} or two numbers. Channels for which tracks are plotted.
%
%    col: {green; red; blue; light grey} or [4x3]-matrix of RGB values.
%         Colours to plot first channel, second channel, delta values, and
%         overlaid plots, respectively.
%
%    filter: {0} or 1. Filter plots by either z-depth, or z-angle < 30deg.
%
%    minLength: {0.25} or number in range 0 to 1. The proportion of a track
%         that needs to be successfully tracked in order to allow plotting.
%
%    plotCoord: {1}, 2 or 3. Coordinate for which trajectories are plotted:
%               1=x, 2=y, 3=z.
%
%    plots: {1} or some subset of the [1:3]. Designates which data to
%           plot:
%               1 = sister pair trajectories
%               2 = interkinetochore distances
%               3 = intrakinetochore distances (including their average)
%
%    subset: {all sisters} or some vector of sister numbers. Subset of
%            sisters for plotting.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2014 C. A. Smith

if nargin<1 || isempty(job)
    error('Need to supply a job.')
end

% set up default options
opts.channels = 1:2;
opts.col = [ 0 ,0.5, 0  ;
             1 , 0 , 0  ;
             0 , 0 , 1 ];
opts.filter = 1;
opts.filterAngle = 30;
opts.filterDepth = 0.1;
opts.filterType = 'zDepth';
opts.minLength = 0.25;
opts.plotCoord = 1;
opts.plots = 1;
opts.subset = [];

% and process the input options
opts = processOptions(opts, varargin{:});

% get initCoord and trackList information
dS1 = job.dataStruct{opts.channels(1)};
dS2 = job.dataStruct{opts.channels(2)};
initCoord1 = dS1.initCoord;
initCoord2 = dS2.initCoord;
trackList1 = dS1.trackList;
trackList2 = dS2.trackList;

% get sisterList information, and ensure sisters are coordinated
sisterList1 = dS1.sisterList;
sisterList2 = dS2.sisterList;
nSisters1 = length(sisterList1);
nSisters2 = length(sisterList2);
if nSisters1 ~= nSisters2
    error('Do not have equal number of sisters. Cannot overlay plots.')
else
    nSisters = nSisters1;
    clear nSisters1 nSisters2
end
% get a list of sister numbers for plotting, if not provided
if isempty(opts.subset)
    opts.subset = 1:nSisters;
elseif max(opts.subset) > nSisters
    opts.subset = 1:nSisters;
    warning('Supplied subset out of range. Will plot all sisters.')
end

% convert filterType to 'none' if no filter required
if ~opts.filter
    opts.filterType = 'none';
end

% find which sisters should not be plotted to avoid unnecessary processing
newSubset = [];
% find maximum length of trajectories
trajLength = job.metadata.nFrames;
for iSis = opts.subset
    
    % calculate for individual scenarios
    failed1 = sum(~isnan(sisterList1(iSis).coords1(:,opts.plotCoord)))/trajLength <= opts.minLength;
    failed2 = sum(~isnan(sisterList1(iSis).coords2(:,opts.plotCoord)))/trajLength <= opts.minLength;
    failed3 = sum(~isnan(sisterList2(iSis).coords1(:,opts.plotCoord)))/trajLength <= opts.minLength;
    failed4 = sum(~isnan(sisterList2(iSis).coords2(:,opts.plotCoord)))/trajLength <= opts.minLength;

    % combine for ultimate failure
    failed = (failed1*failed2 + failed3*failed4)~=0;
    % if not failed, add sister to the list
    if ~failed
        newSubset = [newSubset,iSis];
    end

end
opts.subset = newSubset;

% process plot labels
switch opts.plotCoord
    case 1
        yLabPlots = 'x position';
    case 2
        yLabPlots = 'y position';
    case 3
        yLabPlots = 'z position';
    otherwise
        yLabPlots = '';
end

% plot sister tracks
if any(opts.plots == 1)
   
    % make figure environment
    figure(); clf;
    % get the time information
    t = job.metadata.frameTime(1,:);
    % find the minimum shape of the subplots
    n = length(opts.subset);
    fig_n=ceil(sqrt(n));
    fig_m=ceil(n/fig_n);
    
    % loop over sisters
    for iSis = 1:length(opts.subset)
        % find specific subplot
        subplot(fig_m,fig_n,iSis);
        hold on
        % get trackID and featID information
        trackIDs = sisterList1(1).trackPairs(opts.subset(iSis),1:2);
        % preallocate microscope coords structures for collating initCoords
        micrSisterList1 = sisterList1;
        micrSisterList2 = sisterList2;
        
        featIDs(:,1:2) = [trackList1(trackIDs(1)).featIndx trackList1(trackIDs(2)).featIndx];
        for iTime = 1:trajLength
            % get the featIDs for this kinetochore pair at this timepoint,
            % then transfer initCoord information to the new sisterList
            thisID = featIDs(iTime,:);
            if isnan(thisID(1))
                micrSisterList1(iSis).coords1(iTime,1:3) = NaN;
                micrSisterList2(iSis).coords1(iTime,1:3) = NaN;
            else
                micrSisterList1(iSis).coords1(iTime,1:3) = initCoord1(iTime).allCoord(thisID(1),1:3);
                micrSisterList2(iSis).coords1(iTime,1:3) = initCoord2(iTime).allCoord(thisID(1),1:3);
            end
            if isnan(thisID(2))
                micrSisterList1(iSis).coords2(iTime,1:3) = NaN;
                micrSisterList2(iSis).coords2(iTime,1:3) = NaN;
            else
                micrSisterList1(iSis).coords2(iTime,1:3) = initCoord1(iTime).allCoord(thisID(2),1:3);
                micrSisterList2(iSis).coords2(iTime,1:3) = initCoord2(iTime).allCoord(thisID(2),1:3);
            end
        end
        
        % get sisterList coordinates for the given pair
        pair1 = sisterList1(opts.subset(iSis));
        pair2 = sisterList2(opts.subset(iSis));
        
        switch opts.filterType
        
            case 'zAngle'
                % calculate z-angle
                zRotDiff(1,:) = micrSisterList2(iSis).coords1(:,3)-micrSisterList1(iSis).coords1(:,3);
                xyRotDiff(1,:) = sqrt(sum(micrSisterList2(iSis).coords1(:,1:2)-micrSisterList1(iSis).coords1(:,1:2),2));
                zAngle(1,:) = atan(zRotDiff(1,:)./xyRotDiff(1,:))*180/pi;

                zRotDiff(2,:) = micrSisterList2(iSis).coords2(:,3)-micrSisterList1(iSis).coords2(:,3);
                xyRotDiff(2,:) = sqrt(sum(micrSisterList2(iSis).coords2(:,1:2)-micrSisterList1(iSis).coords2(:,1:2),2));
                zAngle(2,:) = atan(zRotDiff(2,:)./xyRotDiff(2,:))*180/pi;

                % find points being kept after filtering
                filtered(:,1) = +(abs(zAngle(1,:))<opts.filterAngle);
                filtered(:,2) = +(abs(zAngle(2,:))<opts.filterAngle);
                filtered(filtered==0) = NaN;
                
            case 'zDepth'

                % calculate z-depth for filtering if required
                zRotDiff(1,:) = micrSisterList2(iSis).coords1(:,3)-micrSisterList1(iSis).coords1(:,3);

                zRotDiff(2,:) = micrSisterList2(iSis).coords2(:,3)-micrSisterList1(iSis).coords2(:,3);

                % find points being kept after filtering
                filtered(:,1) = +(abs(zRotDiff(1,:))<opts.filterDepth);
                filtered(:,2) = +(abs(zRotDiff(2,:))<opts.filterDepth);
                filtered(filtered==0) = NaN;

            otherwise
                
                % if no filtering required, produce vector of 1s
                filtered = ones(job.metadata.nFrames,2);

        end
        
        % plot the trajectories
        plot(t, pair1.coords1(:,opts.plotCoord).*filtered(:,1),'Color',opts.col(1,:));
        plot(t, pair1.coords2(:,opts.plotCoord).*filtered(:,2),'Color',opts.col(1,:));
        plot(t, pair2.coords1(:,opts.plotCoord).*filtered(:,1),'Color',opts.col(2,:));
        plot(t, pair2.coords2(:,opts.plotCoord).*filtered(:,2),'Color',opts.col(2,:));
        
        % aesthetics per subplot
        sisTit = sprintf('%i',opts.subset(iSis));
        title(sisTit)
        if mod(iSis-1,fig_n)==0
            ylabel(yLabPlots);
        end
        if iSis>(fig_m-1)*fig_n
            xlabel('Time (secs)');
        end
        xlim([min(t) max(t)]);
    end
end

% plot 3D sister separations for both channels
if any(opts.plots == 2)
    figure();
    clf;
    for iSis=1:length(opts.subset)
        subplot(fig_m,fig_n,iSis);
        sep1 = sisterList1(opts.subset(iSis)).distances.*prod(filtered,2);
        sep2 = sisterList2(opts.subset(iSis)).distances.*prod(filtered,2);
        hold on
        plot(t, sep1(:,opts.plotCoord),'Color',opts.col(1,:));
        plot(t, sep2(:,opts.plotCoord),'Color',opts.col(2,:));
        if mod(iSis-1,fig_n)==0
            ylabel('sisterSep');
        end
        if iSis>(fig_m-1)*fig_n
            xlabel('time, s');
        end
        xlim([min(t) max(t)]);
    end
end

% plot 1D and 3D delta
if any(opts.plots == 3)
    figure();
    clf;
    for iSis=1:length(opts.subset)
        
        subplot(fig_m,fig_n,iSis);
        pair1 = sisterList1(opts.subset(iSis));
        pair2 = sisterList2(opts.subset(iSis));
        
        % calculate 1D and 3D delta
        avgDel = abs(sisterList2(opts.subset(iSis)).distances - sisterList1(opts.subset(iSis)).distances).*prod(filtered,2)/2;
        delta1 = sqrt(sum(pair1.coords1(:,1:3) - pair2.coords1(:,1:3),2)).*prod(filtered,2);
        delta2 = sqrt(sum(pair1.coords2(:,1:3) - pair2.coords2(:,1:3),2)).*prod(filtered,2);
        
        hold on
        plot(t, avgDel(:,1),'Color',opts.col(3,:));
        plot(t, delta1(:,1),'Color',opts.col(1,:));
        plot(t, delta2(:,1),'Color',opts.col(2,:));
        if mod(iSis-1,fig_n)==0
            ylabel('delta');
        end
        if iSis>(fig_m-1)*fig_n
            xlabel('time, s');
        end
        xlim([min(t) max(t)]);
    end
end




