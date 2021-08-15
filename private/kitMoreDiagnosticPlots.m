function kitMoreDiagnosticPlots(job,channel,ableToSave)
%This function take a job struct and produces diagnostic plots similar to
%appendix of Jaqaman 2010 to show average distance between sisters and
%average displacement.
%Jonathan U Harrison 2019-02-25
%%%%%%%%%%%%%%%%

if nargin<2
    channel=1;
end
if nargin<3
    ableToSave=true;
elseif ~isfield(job,'output')
    warning('job does not have output to save plots so will not save')
    ableToSave=false;
end
if isfield(job.dataStruct{channel},'sisterList') && ...
        ~isempty(job.dataStruct{channel}.sisterList(1).distances)
    sisterList = job.dataStruct{channel}.sisterList;
    sisterSisterDist = [];
    for i = 1:length(sisterList)
        sisterSisterDist = [sisterSisterDist; sisterList(i).distances(...
            ~isnan(sisterList(i).distances(:,1)),:)];
    end

    figure; histogram(sisterSisterDist(:,1),'Normalization','probability');
    xlabel('Average sister separation ($\mu$m)','Interpreter','Latex');
    ylabel('Normalized frequency','Interpreter','Latex');
    set(gca,'fontsize',20)
    if ableToSave
        savename = sprintf('%sSisterSpearation',kitGenerateOutputFilename(job));
        print(sprintf('%s.eps',savename),'-depsc');
    end
else
    warning('sisterList not found. Unable to complete diagnostic plot');
end
if isfield(job.dataStruct{channel},'trackList')
    %%%%%%%%%%%%
    %FRAME TO FRAME DISPLACEMENTS
    %INCLUDE ALL TRACKS
    trackList = job.dataStruct{channel}.trackList;
    framewiseDistStore = [];
    for i =1:length(trackList)
        framewiseDiff = diff(trackList(i).coords);
        framewiseDist = sqrt(sum(framewiseDiff(:,1:3).^2,2));
        framewiseDistStore = [framewiseDistStore;framewiseDist(~isnan(framewiseDist))];
    end

    figure; histogram(framewiseDistStore,'Normalization','probability');
    xlabel('Frame-to-frame displacement ($\mu$m)','Interpreter','Latex');
    ylabel('Normalized frequency','Interpreter','Latex');
    set(gca,'fontsize',20)
    if ableToSave
        savename = sprintf('%sFramewiseDist',kitGenerateOutputFilename(job));
        print(sprintf('%s.eps',savename),'-depsc');
    end
else
    warning('trackList not found. Unable to complete diagnostic plot');
end


