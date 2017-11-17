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
%    plotAx: {[1,2,3]} or subset of the default. Plot the specified axes,
%    where 1=x, 2=y, 3=z.
%
%    subset: {all tracks} or some vector of tracks. Subset of tracks for plotting.
%
%    minLength: {0.25} or number. Minimum number of tracked frames. Overridden by subset option.
%
% Created by: Jonathan W. Armond 2013
% Edited by:  Chris Smith 10/2013

if nargin<1
    error('Must supply JOB');
end

% Set defaults
opts.channel = 1;
opts.cutoff = 6;
opts.overlay = 1;
opts.plotAx = 1:3;
opts.subset = [];
opts.plotPole = 0;
opts.nLongest = 0;
opts.minLength = 0.25;
% Process options
opts = processOptions(opts, varargin{:});

t = job.metadata.frameTime;
dt = t(1,2)-t(1,1);

dataStruct = job.dataStruct{opts.channel};
trackList = dataStruct.trackList;
nTracks = length(trackList);
maxTime = dataStruct.dataProperties.movieSize(4);

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
n=length(opts.subset);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
clf;
hold on
axName = ['x','y','z'];

for j=1:length(opts.subset)

  i = opts.subset(j);

  x1 = trackList(i).coords(:,1);
  if abs(nanmean(x1))>opts.cutoff && sum(isnan(x1))<ceil(0.5*maxTime);
      poleSub = [poleSub;i];
  end
  t = ((1:length(x1))-1)*dt;

  if opts.overlay == 0
      subplot(fig_m,fig_n,j);
      plot(t,x1);
      title(['Track ' num2str(i)]);
      xlim([0 max(t)])
      xlabel('time, s'); ylabel('x-distance, µm');
      set(gca,'FontSize',8)
  else
      for h=opts.plotAx
          subplot(length(opts.plotAx),1,h)
          hold on
          title(axName(h))
          x1=trackList(i).coords(:,h);
          plot(t,x1);
          xlabel('time, s'); ylabel('x-distance, µm');
          ylim([-12 12])
          set(gca,'FontSize',20)
      end
  end

end

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
        plot(t,x1);
        xlabel('time, s'); ylabel('x-distance, µm');
        set(gca,'FontSize',20)

    end

end
