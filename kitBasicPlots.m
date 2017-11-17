function kitBasicPlots(job, channel)
%KITBASICPLOTS Produces a set of basic plots for a job
%
% SYNOPSIS: kitBasicPlots(job)
%
% INPUT job: Struct containing tracking job setup options.
%
% Copyright (c) 2012 Jonathan W. Armond

% FIXME Call kitPlotSisters, etc..

ds = job.dataStruct{channel};
sisterList = ds.sisterList;
nSisters = length(sisterList);

% Plot sister tracks in X.
figure(1);
t = job.metadata.frameTime(1,:);
n = nSisters;
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
clf;
for i=1:n
    subplot(fig_m,fig_n,i);
    pair = sisterList(i);
    plot(t, pair.coords1(:,1), t, pair.coords2(:,1));
    if mod(i-1,fig_n)==0
        ylabel('X position');
    end
    if i>(fig_m-1)*fig_n
        xlabel('Time (secs)');
    end
    xlim([min(t) max(t)]);
end

