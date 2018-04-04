function kitPlotSisters(job,varargin)
% KITPLOTSISTERS Plot sister tracks as a diagnostic
%
% 'minLength' : set to percentage of track length as threshold for plotting.

if nargin<1
    error('Must supply JOB');
end

% Set defaults
opts.channel = 1;
opts.subset = [];
opts.useTracks = 0;
opts.minLength = 0.75;
% Process options
opts = processOptions(opts, varargin{:});

t = job.metadata.frameTime;
dt = t(1,2)-t(1,1);

dataStruct = job.dataStruct{opts.channel};
if dataStruct.failed || ~isfield(dataStruct,'sisterList')
  fprintf('\nSpot finding failed for this movie.\n\n');
  return
end
sisterList = dataStruct.sisterList;
if isempty(sisterList(1).trackPairs)
  fprintf('\nNo sisters found in this movie.\n\n');
  return
end

trackList = dataStruct.trackList;
trackPairs = sisterList(1).trackPairs;
if opts.minLength > 0 && isempty(opts.subset)
  coords = horzcat(sisterList.coords1);
  coords = coords(:,1:6:end); % X coordinate.
  nancount = sum(isnan(coords),1);
  opts.subset = find(nancount <= job.metadata.nFrames*(1-opts.minLength));
end
if isempty(opts.subset)
  opts.subset = 1:length(sisterList);
end
sisterList = sisterList(opts.subset);
nSisters = length(sisterList);

figure;
n=nSisters;
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
clf;
for i=1:nSisters
    subplot(fig_m,fig_n,i);
    pair = trackPairs(i,1:2);
    if opts.useTracks
      x1 = trackList(pair(1)).coords(:,1);
      x2 = trackList(pair(2)).coords(:,1);
    else
      x1=sisterList(i).coords1(:,1);
      x2=sisterList(i).coords2(:,1);
    end
    t=((1:length(x1))-1)*dt;
    plot(t,x1,t,x2);
    title(sprintf('sister %d tracks %d,%d',opts.subset(i),pair(1),pair(2)));
    xlim([t(1) t(end)]);
    xlabel('t');
    ylabel('x');
end
