function maxSpots = kitCheckTracks(job,saveOutput)
% As a diagnostic tool it is helpful to investigate the number of tracks
% covering three consecutive frames
%
% saveOutput - boolean to indicate whether to save plots. If false, will close all open figures. 
%
%
%Jonathan U Harrison 2019-02-21
%%%%%%%%%%%%%%%%%
if nargin < 2
    saveOutput = 1;
end

if isfield(job,'dataStruct') && isfield(job.dataStruct{1},'tracks')
    %compute how many tracks in a given frame
    nTracks = size(job.dataStruct{1}.tracks,1);
    nFrames = job.metadata.nFrames;
    nStages = 2;
    nSpots = zeros(nFrames,nStages);
    t = job.metadata.frameTime(1,:);
    if isfield(job,'output')
        savename =  sprintf('%s',job.output);
        savename(regexp(savename,'[.]'))=[];
        ableToSave = true;
    else
        warning('unable to find job.output field so unable to save');
        ableToSave = false;
    end
    
    
    doesTrackExist = zeros(nFrames,nTracks);
    figure; hold on; xlabel('Time (s)'); ylabel('Track ID');
    title('Barcode plot showing birth and death of tracklets');
    set(gca,'fontsize',20);
    longEnoughCutoff = job.options.minSisterTrackOverlap; 
    for j = 1:nTracks
        %put ones where a track exists
        startTime = job.dataStruct{1}.tracks(j).seqOfEvents(1,1);
        endTime = job.dataStruct{1}.tracks(j).seqOfEvents(2,1);
        doesTrackExist(startTime:endTime,j) = ones(endTime-startTime+1,1);
        if endTime - startTime>longEnoughCutoff
            barcodeColour = 'r';
        else
            barcodeColour = 'k';
        end
        plot(t(startTime:endTime),j*ones(size(startTime:endTime)),...
            [barcodeColour,'-']);
    end
    if ableToSave & saveOutput
        print(sprintf('%sBarcode.eps',savename),'-depsc');
    end
    nSpots(:,1) = sum(doesTrackExist,2);
    maxSpots = max(nSpots(:,1)); %this gives the maximum number of spots that exist at one time (allowing spots that exist for any length of time); potentially using for assessing the number of chromosomes present in each cell

    doesTrackEndHere = zeros(nFrames,nTracks);
    for j = 1:nTracks
        %count how many tracks end at the given frame
        endTime = job.dataStruct{1}.tracks(j).seqOfEvents(2,1);
        doesTrackEndHere(endTime,j) = 1;
    end
    nSpots(:,2) = sum(doesTrackEndHere,2);
    
    %%%%%%%%%%%%%%
    % and now the plotting
    coloursByStage = [230,97,1;
        253,184,99;
        178,171,210;
        94,60,153]/256; %from http://colorbrewer2.org/
    symbolByPlane = 'x';
    stagesLegend = {'track','tracksEnding'};
    figure;
    hold all;
    for i=1:nStages
        scatter(t,nSpots(:,i),100,...
            repmat(coloursByStage(i,:),nFrames,1), symbolByPlane);
    end
    legend(stagesLegend,'Location','southwest');
    xlabel('Time (s)');
    ylabel('Number of spots');
    title('Performance of each stage of KiT in tracking kinetochores');
    grid on;
    box on;
    set(gca,'fontsize',20);
    if ableToSave & saveOutput
        print(sprintf('%sNumSpotsTrackingDiagnosticPlot.eps',savename),'-depsc');
    end
else
    warning('tracks field not found. Unable to perform checks')
    maxSpots = NaN;
end
