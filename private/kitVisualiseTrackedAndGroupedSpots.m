function kitVisualiseTrackedAndGroupedSpots(job,movie,channel,frameToShow)
% Visualise detected and tracked spots before and after grouping the
% sister kinetochore pairs
%
% job - struct, containing job info
% movie - 4-D double, full movie already read in
% channel - int, channel to visualise in
% frameToShow - int, pick a frame to plot sisters on
%
%Jonathan U Harrison 2019-02-22
%%%%%%%%%%%%%%%%
if isfield(job,'metadata') && frameToShow>job.metadata.nFrames
    error('Are you sure that frame exists? Frame requested=%d, but nFrames=%d',...
        frameToShow,job.metadata.nFrames);
end
if isfield(job.dataStruct{channel},'tracks')
    %visualise tracked spots
    totalTrackedSpots = size(job.dataStruct{channel}.tracks,1);
    figure;
    subplot(1,2,1);
    imshow(max(movie(:,:,:,frameToShow),[],3));
    hold on;
    if ~isfield(job.dataStruct{channel},'trackList')
        %if tracks present but not processed, then extract them
        fprintf('Extracting tracks...\n');
        job = kitExtractTracks(job, channel);
        fprintf('Done\n');
    end
    for j=1:totalTrackedSpots
        %loop over tracks and plot at relevant frame if its there
        coords = job.dataStruct{channel}.trackList(j).coords(frameToShow,1:3);
        coords = kitCoordsToImageCoords(job, channel, coords,frameToShow);
        plot(coords(:,1),coords(:,2),'mx');
    end
    pause(0.2);
else
    warning('tracks field not found. Will not perform diagnostic plot');
end

if isfield(job.dataStruct{channel},'sisterList') && ...
        ~isempty(job.dataStruct{channel}.sisterList(1).coords1)
    %visualise grouped sisters
    totalNumSisters = length(job.dataStruct{channel}.sisterList);
    subplot(1,2,2);
    imshow(max(movie(:,:,:,frameToShow),[],3));
    hold on;
    for j=1:totalNumSisters
        coords1=kitCoordsToImageCoords(job, channel, ...
            job.dataStruct{channel}.sisterList(j).coords1(frameToShow,:),...
            frameToShow);
        coords2=kitCoordsToImageCoords(job, channel, ...
            job.dataStruct{channel}.sisterList(j).coords2(frameToShow,:),...
            frameToShow);
        plot(coords1(:,1),coords1(:,2),'co');
        plot(coords2(:,1),coords2(:,2),'rx');
    end
    %%%%%%%%%%%
else
    warning('sisterList field not found. Will not perform diagnostic plot');
end
