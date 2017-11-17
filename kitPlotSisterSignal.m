function kitPlotSisterSignal(job,sigchannel,varargin)
% KITPLOTSISTERSIGNAL Plot sister tracks and signal
%
%    KITPLOTSISTERSIGNAL(JOB,SIGCHANNEL,...) Plots sisters against signal
%    intensity in channel SIGCHANNEL. Supply options as string/value pairs.
%
%    Options, defaults in {}:-
%
%    channel: {1} or number. Tracking channel for which tracks are plotted.
%
%    minLength: {0.75}. set to percentage of track length as threshold for plotting.
%
%    sel: {[]}. subset of tracks.
%
%    directionality: {0}. Set to 1 to plot AP/P sections.
%
% Copyright(c) Jonathan W. Armond 2015

if nargin<2
    error('Must supply JOB and SIGCHANNEL');
end

% Set defaults
opts.channel = 1;
opts.sel = [];
opts.minLength = 0.75;
opts.directionality = 0;

% Process options
opts = processOptions(opts, varargin{:});

ds = job.dataStruct{opts.channel};
sisterList = ds.sisterList;
trackPairs = sisterList(1).trackPairs;
trackList = ds.trackList;
trackInt = ds.trackInt;
t = job.metadata.frameTime(1,:)';

if opts.minLength > 0
  coords = horzcat(sisterList.coords1);
  coords = coords(:,1:6:end); % X coordinate.
  nancount = sum(isnan(coords),1);
  opts.sel = find(nancount < job.metadata.nFrames*(1-opts.minLength));
elseif ~isempty(opts.sel)
  opts.sel = 1:length(sisterList);
end
sisterList = sisterList(opts.sel);
trackPairs = trackPairs(opts.sel,:);
nSisters = length(sisterList);


figure;
n=nSisters;
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
clf;
for i=1:nSisters
  subplot(fig_m,fig_n,i);
  pair = trackPairs(i,1:2);

  % Plot sister tracks.
  plotyy([t t],[sisterList(i).coords1(:,1) sisterList(i).coords2(:,1)],...
         [t t],[trackInt(pair(1)).intensity(:,sigchannel) trackInt(pair(2)).intensity(:,sigchannel)])
  if opts.directionality
    % Plot directionality.
    hold on;
    sc=sisterList(i).coords1(:,1);
    sc(trackList(pair(1)).direction~=+1)=nan;
    plot(t,sc,'g-','linewidth',2);
    sc=sisterList(i).coords2(:,1);
    sc(trackList(pair(2)).direction~=+1)=nan;
    plot(t,sc,'g-','linewidth',2);
    sc=sisterList(i).coords1(:,1);
    sc(trackList(pair(1)).direction~=-1)=nan;
    plot(t,sc,'m-','linewidth',2);
    sc=sisterList(i).coords2(:,1);
    sc(trackList(pair(2)).direction~=-1)=nan;
    plot(t,sc,'m-','linewidth',2);
    hold off;
  end

  xlim([t(1) t(end)]);
  title(sprintf('sister %d tracks %d,%d',i,pair(1),pair(2)));
end
suptitle('KT position blue/red; fluoresence gold/magenta');

% Return data.
%x=[t' sisterList(idx).coords1(:,1) sisterList(idx).coords2(:,1) ...
%   trackInt(pair(1)).intensity_max(:,sigchannel) trackInt(pair(2)).intensity_max(:,sigchannel)];
