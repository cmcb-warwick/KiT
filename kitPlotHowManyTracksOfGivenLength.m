function barcodeData = kitPlotHowManyTracksOfGivenLength(job, channel)
% Given a job, compute barcode data in the form of a
% vector of the duration of tracks. Using this barcode, plot something like
% 1-cdf of length of tracklets to assess quality of tracks
%
%
% Jonathan U Harrison 2019-03-14
%%%%%%%%%%%%%%%%%%

if nargin < 2
    channel=1;
end

if isfield(job,'dataStruct') && isfield(job.dataStruct{channel},'tracks')
fprintf('Loading barcode data from jobset for analysis\n');
%load in barcode data from jobset rather than direct input
nTracks = size(job.dataStruct{channel}.tracks,1);
nFrames = job.metadata.nFrames;
barcodeData = zeros(nTracks,1);
for j = 1:nTracks
    %put ones where a track exists
    startTime = job.dataStruct{channel}.tracks(j).seqOfEvents(1,1);
    endTime = job.dataStruct{channel}.tracks(j).seqOfEvents(2,1);
    barcodeData(j) = (endTime - startTime +1)/nFrames;
end
[F,x,Flo,Fup] = ecdf(barcodeData,'bounds','on'); %requires stats toolbox
x(end)=[]; F(end)=[]; Flo(end)=[]; Fup(end)=[]; %looks better without final pt

plot_reverse_ecdf(nTracks,x,F,Flo,Fup,0);
TrackedSpotsInTracksOfLength(barcodeData,nFrames,0);

%%%%%%%%%%%%%%%%%
%now do the same for grouped sister pairs
nSisters = length(job.dataStruct{channel}.sisterList);
sisterBarcodeData = zeros(nSisters,1);
for i =1:nSisters
   sisterBarcodeData(i) = sum(~(isnan(job.dataStruct{channel}.sisterList(i).coords1(:,1)) & ...
        isnan(job.dataStruct{channel}.sisterList(i).coords2(:,1))))/nFrames;
    %allow either sister 1 or sister 2 to have valid data
end

[F,x,Flo,Fup] = ecdf(sisterBarcodeData,'bounds','on'); %requires stats toolbox
x(end)=[]; F(end)=[]; Flo(end)=[]; Fup(end)=[]; %looks better without final pt
plot_reverse_ecdf(nSisters,x,F,Flo,Fup,1);
TrackedSpotsInTracksOfLength(sisterBarcodeData,nFrames,1);
end
end

function plot_reverse_ecdf(nTracks,x,F,Flo,Fup,useSisters)
nTotal = 46*(2-useSisters);
figure;
plot(x, (ones(size(x)) - F)*nTracks,'k-','linewidth',3);
hold on;
plot(x, (ones(size(x)) - Flo)*nTracks,'k--','linewidth',3);
plot(x, (ones(size(x)) - Fup)*nTracks,'k--','linewidth',3);
plot(x, ones(size(x))*nTotal, 'g.-','linewidth',3);
xlabel('Proportion of full track length')
if ~useSisters
    ylabel('# of tracks of given length');
    title('Number of tracks of a given length');
else %change the labels
    ylabel('# of sisters of given length');
    title('Number of sisters of a given length');
end    
set(gca,'fontsize',20);
ylim([0,100*(2-useSisters)]);
end

function TrackedSpotsInTracksOfLength(barcodeData,nFrames,useSisters)
nTotal = 46*(2-useSisters);
nl = 200;
LArray = linspace(0,1,nl);
S = zeros(nl,1);
for i=1:nl
    S(i) = sum((barcodeData>=LArray(i)).*barcodeData);
end
figure;
plot(LArray,S,'k-','linewidth',3);
hold on;
%plot(x, (ones(size(x)) - Flo)*nTracks,'k--','linewidth',3);
%plot(x, (ones(size(x)) - Fup)*nTracks,'k--','linewidth',3);
plot(LArray, ones(size(LArray))*nTotal, 'g.-','linewidth',3);
xlabel('Proportion of full track length')
if ~useSisters
    ylabel('# tracked spots in tracks of given length');
    title('Number of tracks of a given length');
else %change the labels
    ylabel('# tracked sisters in sister pairs of given length');
    title('Number of sisters of a given length');
end    
set(gca,'fontsize',20);
%ylim([0,100*(2-useSisters)]);
end
