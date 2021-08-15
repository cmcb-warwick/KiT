function kitPlotTracks(job,varargin)
% KITPLOTTRACKS Plot tracks as a diagnostic
%
%    KITPLOTTRACKS(JOB,...) Plots all tracks by default overlaid on
%    separate graphs for each coordinate. Supply options as string/value
%    pairs following JOB.
%
%    Options, defaults in {}:-
%
%    channel: {1} or number. Channel for which tracks are plotted.
%
%    cutoff: 0 or {1}. Define distance, in microns, away from the metaphase
%    plate beyond which a track is a candidate for that of a spindle pole.
%
%    overlay: 0 or {1}. Overlays tracks for each coordinate.
%
%    plotAx: {1} or subset of the default. Plot the specified axes,
%    where 1=x, 2=y, 3=z.
%
%    subset: {all tracks} or some vector of tracks. Subset of tracks for plotting.
%
%    faintSubset: {no tracks} Plot this subset of trajectories faint in the
%    background in grey
%
%    usePairs: {0} or 1. Plot trajectories of sister pairs or individual
%    sisters
%
%    identifyLazyKTs: {0} or 1. Overwrites other subset options. Will
%    determine the lazy or lagging KTs and plot these with others in faint.
%
%    minLength: {0.25} or number. Minimum number of tracked frames. Overridden by subset option.
%
%    dt: {[]} or number. Time step between frames. Otherwise determined from metadata. 
%    Only needed when there is an error in the metadata. Rescales time on x axis if set.
%     
%    overWriteColors: {[]} or array eg. [0.4660    0.6740    0.1880]. Uses this for plotting 
%    colors instead of the otherwise default custom colours.
%
% Created by: Jonathan W. Armond 2013
% Edited by:  Chris Smith 10/2013

if nargin<1
    error('Must supply JOB');
end

% Set defaults
opts.channel = 1;
opts.cutoff = 1000;
opts.overlay = 1;
opts.plotAx = 1;
opts.subset = [];
opts.faintSubset = [];
opts.plotPole = 0;
opts.nLongest = 0;
opts.minLength = 0.25;
opts.usePairs = 0;
opts.identifyLazyKTs = 0;
opts.movie = [];
opts.makeMovie = 1;
opts.dt = [];
opts.overWriteColors = [];
% Process options
opts = processOptions(opts, varargin{:});

t = job.metadata.frameTime;
if isempty(opts.dt)
dt = t(1,2)-t(1,1);
else 
dt = opts.dt;
t = t*dt/2.0; %rescale to fix incorrect metadata
end

dataStruct = job.dataStruct{opts.channel};
trackList = dataStruct.trackList;
nTracks = length(trackList);
maxTime = dataStruct.dataProperties.movieSize(4);

if isempty(dataStruct.sisterList(1).trackPairs)
    fprintf('\nNo sisters found in this movie.\n\n');
if opts.usePairs
    return
else
    fprintf('Atempting to continue with unpaired tracks\n\n');
end
end

if isempty(opts.subset)
    if isempty(opts.minLength)
        opts.subset = 1:nTracks;
    else
        coords = horzcat(trackList.coords);
        coords = coords(:,1:6:end); % X coordinate.
        nancount = sum(isnan(coords),1);
        opts.subset = find(nancount < job.metadata.nFrames*(1-opts.minLength));
    end
end

if opts.nLongest > 0
    n = zeros(length(opts.subset),1);
    for j=1:length(opts.subset)
        n(j) = sum(~isnan(trackList(opts.subset(j)).coords(:,1)));
    end
    [~,idx] = sort(n,'descend');
    opts.subset = opts.subset(idx(1:opts.nLongest));
end

poleSub = [];

figure;
axes('NextPlot', 'add')
clf;
n=length(opts.subset);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
hold on
axName = ['x','y','z'];

if opts.identifyLazyKTs
    lazyKTs = kitIdentifyLazyKTs(job,opts.channel);
    subset = unique(lazyKTs);
    faintSubset = 1:length(job.dataStruct{opts.channel}.trackList);
else
    subset = opts.subset;
    faintSubset = opts.faintSubset;
    if opts.usePairs
        trackPairs = dataStruct.sisterList(1).trackPairs(:,1:2);
        subset = trackPairs(subset,1:2);
        faintSubset = trackPairs(faintSubset,1:2); faintSubset = faintSubset(:)';
    else
        subset = subset';
    end
end

if ~isempty(faintSubset)
    for k=faintSubset
        for h=opts.plotAx
            x1=trackList(k).coords(:,h);
            t = ((1:length(x1))-1)*dt;
            subplot(length(opts.plotAx),1,h)
            title(axName(h))
            plot(t,x1,'color',[0 0 0 0.1],'linewidth',1);
            xlim([0 max(t)])
            xlabel('Time, s'); ylabel('Position, µm');
            ylim([-12 12]);
            set(gca,'FontSize',20);
            hold on
        end
    end
else
    % for some reason the interactive plot needs subplots already made,
    % so do this even if not using the faint subset
    for h=opts.plotAx
        subplot(length(opts.plotAx),1,h)
        title(axName(h))
        xlabel('Time, s'); ylabel('Position, µm');
        set(gca,'FontSize',20);
        hold on
    end
end

transparentplt = [];
for ii =1:(opts.usePairs+1)
    [poleSub,plt] = plot_subset_highlighted_tracks(subset(:,ii)',trackList,opts,...
        dt,maxTime,poleSub,fig_n,fig_m,axName,job);
end
transparentplt = [transparentplt,plt];

if opts.plotPole && ~isempty(poleSub)
    figure;
    n=length(poleSub);
    fig_n=ceil(sqrt(n));
    fig_m=ceil(n/fig_n);
    clf;
    hold on
    
    for j=1:length(poleSub)
        
        i=poleSub(j);
        x1 = trackList(i).coords(:,1);
        
        subplot(fig_m,fig_n,j);
        title(['Track ' num2str(i)]);
        pause(1)
        plot(t,x1,'linewidth',3);
        xlim([0 max(t)])
        xlabel('time, s'); ylabel('x-position, µm');
        set(gca,'FontSize',20)
        
    end
    
end
end

function [poleSub, transparentplt] = plot_subset_highlighted_tracks(subset,trackList,opts,...
    dt,maxTime,poleSub,fig_n,fig_m,axName,job)
ColorOdrCustom = [0 0 1;...
    0 1 0;...
    1 0 0;...
    0 1 1;...
    1 0 1;...
    1 0.69 0.39;...
    0.6 0.2 0;...
    0 0.75 0.75;...
    0.22 0.44 0.34;...
    0.32 0.19 0.19]; %10x3 RGB array See https://uk.mathworks.com/matlabcentral/answers/133676-change-automatically-colors-and-$
if ~isempty(opts.overWriteColors)
ColorOdrCustom = repmat(opts.overWriteColors,size(ColorOdrCustom,1),1);
end

for j=1:length(subset)
    
    i = subset(j);
    x1 = trackList(i).coords(:,1);
    t = ((1:length(x1))-1)*dt;
    if abs(nanmean(x1))>opts.cutoff && sum(isnan(x1))<ceil(0.5*maxTime);
        poleSub = [poleSub;i];
    end
    
    if opts.overlay == 0
        subplot(fig_m,fig_n,j);
        plot(t,x1,'color',ColorOdrCustom(mod(j,10)+1,:),'linewidth',3);
        title(['Track ' num2str(i)]);
        xlim([0 max(t)])
        xlabel('Time, s'); ylabel('x-position, µm');
        set(gca,'FontSize',20)
    else
        for h=opts.plotAx
            subplot(length(opts.plotAx),1,h)
            hold on
            title(axName(h))
            x1=trackList(i).coords(:,h);
            transparentplt(j) = plot(t,x1,'color',ColorOdrCustom(mod(j,10)+1,:),'linewidth',3);
            xlim([0 max(t)])
            xlabel('Time, s'); ylabel('Position, µm');
            ylim([-12 12])
            set(gca,'FontSize',20)
        end
    end
    
end
set(transparentplt, 'ButtonDownFcn', {@LineSelected, transparentplt, trackList, job, opts})
end
function LineSelected(ObjectH, EventData, H, trackList, job, opts)
% This function provides interactive functionality to slecet a particular
% track and visualise this on the movie.
% Basic function for selecting lines based on https://uk.mathworks.com/matlabcentral/answers/83231-how-to-use-the-mouse-to-select-and-identify-a-line-on-a-plot?s_tid=answers_rc1-2_p2_Topic
% Jonathan U. Harrison 2019
%%%%%%%%%%%%%%%%%%

set(ObjectH, 'LineWidth', 5);
set(H(H ~= ObjectH), 'LineWidth', 3);
for i =1:length(trackList)
    if nansum(abs(ObjectH.YData' - trackList(i).coords(:,1))) < 10^(-12) %test for equality with small tolerance
        fprintf('You have selected the %dth track. Plotting this on movie... \n',i);
        if isempty(opts.movie)
            fprintf('No movie provided so loading movie ...\n')
            if isfield(job,'ROI') && length(job.ROI)>1
                warning('More than one stack in job but using only first one\n');
                pathToMovie = fullfile(job.movieDirectory,job.ROI(1).movie);
            else
                pathToMovie = fullfile(job.movieDirectory,job.ROI.movie);
            end
            if isfield(job,'metadata')
                if iscell(job.metadata)
                    job.metadata = job.metadata{job.index};
                end
                [job.metadata, reader] = kitOpenMovie(pathToMovie,'valid',job.metadata);
            else
                [job.metadata, reader] = kitOpenMovie(pathToMovie);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %load movie
            movie = kitReadWholeMovie(reader,job.metadata,opts.channel,job.ROI(1).crop,0,1);
            fprintf('Done loading\n');
        else
            movie = opts.movie;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        trackIDs = job.dataStruct{opts.channel}.sisterList(1).trackPairs(:,1:2);
        refChan = job.options.coordSystemChannel;
        chrShift = job.options.chrShift.result{refChan,opts.channel}(1:3);
        pixelSize = job.metadata.pixelSize;
        makeMovie = opts.makeMovie;
            hh = figure;
            if ~verLessThan('matlab','9.5') %this property only added for matlab R2018a onwards
                hh.WindowState = 'fullscreen'; %make figure full screen
            end
        if makeMovie
            outname = kitGenerateOutputFilename(job);
            vidfile = VideoWriter(sprintf('%sSelectedTrackMovie_Track%d.mp4',outname,i),'MPEG-4');
            open(vidfile);
        end
        if any(i==trackIDs(:))
            fprintf('Selected track is paired. Plotting both\n' );
            for tP = 1:job.metadata.nFrames %loop over time
                coords = nan(2,3); %coords x sister
                [sisID,~] = ind2sub(size(trackIDs),find(i==trackIDs(:))); %allows to find corresponding sister pair
                
                % accumulate track information by channel and sister
                for iSis = 1:2
                    tk = trackIDs(sisID,iSis);
                    track = job.dataStruct{opts.channel}.tracks(tk);
                    
                    startTime = track.seqOfEvents(1,1);
                    endTime   = track.seqOfEvents(2,1);
                    if tP < startTime || tP > endTime
                        coords(iSis,:) = nan(1,3);
                    else
                        coords(iSis,:) = ...
                            track.tracksCoordAmpCG(8*(tP-(startTime-1))-7:8*(tP-(startTime-1))-5);
                        coords(iSis,:) = coords(iSis,:) + chrShift;
                        coords(iSis,:) = coords(iSis,:)./pixelSize;
                    end
                end
                imshow(max(movie(:,:,:,tP),[],3),[]); hold on; plot(coords(:,1),coords(:,2),'rx','markerSize',12);
                if makeMovie
                    drawnow;
                    F(tP) = getframe(hh);
                    writeVideo(vidfile,F(tP));
                end
            end
        else
            fprintf('Selected track is not paired. Plotting individually\n');
            for tP = 1:job.metadata.nFrames %loop over time
                coords = nan(1,3);
                tk = i;
                track = job.dataStruct{opts.channel}.tracks(tk);
                startTime = track.seqOfEvents(1,1);
                endTime   = track.seqOfEvents(2,1);
                
                if tP < startTime || tP > endTime
                    coords = nan(1,3);
                else
                    coords(:) = ...
                        track.tracksCoordAmpCG(8*(tP-(startTime-1))-7:8*(tP-(startTime-1))-5);
                    coords = coords + chrShift;
                    coords = coords./pixelSize;
                end
                imshow(max(movie(:,:,:,tP),[],3),[]); hold on; plot(coords(1),coords(2),'rx','markerSize',12);
                if makeMovie
                    drawnow;
                    F(tP) = getframe(hh);
                    writeVideo(vidfile,F(tP));
                end
            end
            if makeMovie
                close(vidfile)
            end
        end
    end
end
end
