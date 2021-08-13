function kitPlateCheck(job,varargin)
% KITPLATECHECK Visualize plate in yz-plane for coord-sys verification
%
%    KITPLATECHECK(JOB,...) Visualize plate in xy-plane for
%    coordinate-system verification. JOB may be either a job struct, loaded via
%    kitLoadJob, or a sisterList.
%
%    'channel': Optional, specify tracked channel to display.
%
%    'colors': Optional, set to 1 to show multiple colors to distinguish
%    tracks. Otherwise uses colors to distinguish sisters
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
%    Jonathan U. Harrison 2019-06-25
%%%%%%%%%%%%%%%%%%%%%

opts.channel = 1;
opts.subset = [];
opts.colors = 0;
opts.faintSubset = [];
opts.usePairs = 0;
opts.identifyLazyKTs = 0;
opts.minLength = 0.25;
opts = processOptions(opts,varargin{:});

if isfield(job,'dataStruct')
    % Is a job struct.
    sisterList = job.dataStruct{opts.channel}.sisterList;
    dataStruct = job.dataStruct{opts.channel};
    trackList = dataStruct.trackList;
    nTracks = length(trackList);
else
    % Assume a sisterList.
    error('this will no longer work with a sisterList');
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

ColorOdrCustom = [0 0 1;...
    0 1 0;...
    1 0 0;...
    0 1 1;...
    1 0 1;...
    1 0.69 0.39;...
    0.6 0.2 0;...
    0 0.75 0.75;...
    0.22 0.44 0.34;...
    0.32 0.19 0.19]; %10x3 RGB array See https://uk.mathworks.com/matlabcentral/answers/133676-change-automatically-colors-and-line-style

coords1 = horzcat(sisterList.coords1);
coords2 = horzcat(sisterList.coords2);

if opts.identifyLazyKTs
    lazyKTs = kitIdentifyLazyKTs(job,opts.channel);
    subset = unique(lazyKTs);
    faintSubset = 1:length(job.dataStruct{opts.channel}.trackList);
else
    subset = opts.subset;
    faintSubset = opts.faintSubset;
    if opts.usePairs
        trackPairs = dataStruct.sisterList(1).trackPairs(:,1:2);
        subset = trackPairs(subset,1:2); subset = subset(:)';
        faintSubset = trackPairs(faintSubset,1:2); faintSubset = faintSubset(:)';
    end
end
figure;
for k=faintSubset
    y1=trackList(k).coords(:,2);
    z1=trackList(k).coords(:,3);
    plot(y1,z1,'color',[0 0 0 0.1],'linewidth',1);
    xlabel('Y Position (um)');
    ylabel('Z Position (um)');
    ylim([-12 12]);
    set(gca,'FontSize',20);
    hold on
end

if opts.usePairs
    for j=1:length(subset)/2
        i = subset(j);
        iPair = subset(j+length(subset)/2);
        y1 = trackList(i).coords(:,2);
        z1 = trackList(i).coords(:,3);
        hold on %note both in pair should be the same color
        plot(y1,z1,'color',ColorOdrCustom(mod(j,10)+1,:),'linewidth',3);
        y1 = trackList(iPair).coords(:,2);
        z1 = trackList(iPair).coords(:,3);
        %make pairs similar but slightly different colours
        plot(y1,z1,'color',[max(min(ColorOdrCustom(mod(j,10)+1,:)+0.5*randn(1,3),1),0),0.5],'linewidth',3);
        xlabel('Y Position (um)');
        ylabel('Z Position (um)');
        ylim([-12 12])
        set(gca,'FontSize',20)
    end
else
    for j=1:length(subset)
        
        i = subset(j);
        
        y1 = trackList(i).coords(:,2);
        z1 = trackList(i).coords(:,3);
        hold on
        plot(y1,z1,'color',ColorOdrCustom(mod(j,10)+1,:),'linewidth',3);
        xlabel('Y Position (um)');
        ylabel('Z Position (um)');
        ylim([-12 12])
        set(gca,'FontSize',20)
    end
end
%figure;
%if opts.colors
%  plot(coords1(:,2:6:end),coords1(:,3:6:end),...
%       coords2(:,2:6:end),coords2(:,3:6:end));
%else
%  plot(coords1(:,2:6:end),coords1(:,3:6:end),'b-',...
%       coords2(:,2:6:end),coords2(:,3:6:end),'g-');
%end
%xlabel('y');
%ylabel('z');
