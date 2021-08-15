function kitShowGroupedPairs(job,channel,opts)
%%name could be improved
%%
%% plot crosses on an image corresponding only to the grouped sister pairs
%%
%% Jonathan U Harrison 2019-11-07
%%%%%%%%%%%%%%
if nargin<3
    opts.transpose=0; %whether to flip the image and coordinates
    opts.makeMovie=1; %save the output as an mp4 file
    opts.joinSisters=1;
end
if nargin<2
    channel = 1;
end
if ~isfield(job,'dataStruct')
    error('No dataStruct found in job. Did tracking work properly?');
end
%read in movie
if isfield(job,'metadata')
    if iscell(job.metadata)
        job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
end
nTimePoints = job.metadata.nFrames;

if opts.makeMovie
    outname = kitGenerateOutputFilename(job);
    vidfile = VideoWriter(sprintf('%sPairedSistersMovie.mp4',outname),'MPEG-4');
    open(vidfile);
end
h=figure;
if ~verLessThan('matlab','9.5') %this property only added for matlab R2018a onwards
    h.WindowState = 'fullscreen'; %make figure full screen
end

for timePoint=1:nTimePoints
    img = kitReadImageStack(reader,job.metadata,timePoint,channel,job.ROI.crop([2,1,4,3]),0);
    if opts.transpose
        img = img';
    end
    imshow(max(img,[],3),[]);
    hold on;
    for sisPair=1:length(job.dataStruct{channel}.sisterList)
        coords = getTrackInfo(job,channel,sisPair,timePoint);
        if opts.transpose
            for iSis=1:2
            plot(coords(iSis,2),coords(iSis,1),...
                'Color','r','Marker','x','MarkerSize',15);
            hold on;
            end
            %add a line between sisters
            plot(coords(:,2),coords(:,1),...
                'Color','r');
            hold on;
        else
            for iSis=1:2
            plot(coords(iSis,1),coords(iSis,2),...
                'Color','r','Marker','x');
            end
            if opts.joinSisters
            %add a line between sisters
            plot(coords(:,1),coords(:,2),...
                'Color','r');
            hold on;
            end
        end
    end
    title(sprintf("Frame: %d",timePoint));
    set(gca,'fontsize',20);
    pause(0.05);
    if opts.makeMovie
        drawnow;
        F(timePoint) = getframe(h);
        writeVideo(vidfile,F(timePoint));
    end
end
if opts.makeMovie
    close(vidfile)
end
%%%%%%%%%%%%%%%%%
% h = figure;
% if ~verLessThan('matlab','9.5') %this property only added for matlab R2018a onwards
%     h.WindowState = 'fullscreen'; %make figure full screen
% end
% outname = kitGenerateOutputFilename(job);
% vidfile = VideoWriter(sprintf('%sPairedSistersMovie.mp4',outname),'MPEG-4');
% open(vidfile);
% for i = 1:nTimePoints
%     showPlaneFit(job,reader,channel,i,planeFit(i).planeOrigin,...
%         planeFit(i).planeVectors,initCoord(i).allCoordPix);
%     if opts.debug.makePlaneFitMovie
%         drawnow;
%         F(i) = getframe(h);
%         writeVideo(vidfile,F(i));
%     end
% end
% if opts.debug.makePlaneFitMovie
%     close(vidfile)
% end

end % function kitFitPlane

%% LOCAL FUNCTIONS

function coords = getTrackInfo(job,channel,sisPair,timePoint)

sisterList = job.dataStruct{channel}.sisterList;
% get track information
trackIDs = sisterList(1).trackPairs(sisPair,1:2);
%%%%%%%%%%%%

% get pixel resolution
pixelSize = job.metadata.pixelSize;
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,channel}(1:3);

coords = nan(2,3); %coords x sister

% accumulate track information by channel and sister
for iSis = 1:2
    tk = trackIDs(iSis);
    track = job.dataStruct{channel}.tracks(tk);
    
    startTime = track.seqOfEvents(1,1);
    endTime   = track.seqOfEvents(2,1);
    if timePoint < startTime || timePoint > endTime
        coords(iSis,:) = nan(1,3);
    else
        coords(iSis,:) = ...
            track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
        coords(iSis,:) = coords(iSis,:) + chrShift;
        coords(iSis,:) = coords(iSis,:)./pixelSize;
    end
end
end
% 
% function showPlaneFit(job,reader,channel,frameNum,origin,eVecs,allCoordPix)
% % Show debugging visualization.
% 
% img = kitReadImageStack(reader, job.metadata, frameNum, channel, job.ROI.crop);
% % Max project image.
% img = max(img,[],3);
% 
% % Draw image.
% imshow(img,[]);
% title(sprintf('Plane fit frame %d/%d (x red, y cyan, z yellow)',...
%     frameNum,job.metadata.nFrames));
% hold on;
% 
% % Draw axes and origin. FIXME only looking at 2D fit
% origin = origin ./ job.metadata.pixelSize;
% origin = origin(1:2)';
% plot(origin(1),origin(2),'go');
% if ~isempty(eVecs)
%     axisLen = 30;
%     eVecs = eVecs(1:2,1:3);
%     xAxis = [origin origin+axisLen*eVecs(:,1)];
%     yAxis = [origin origin+axisLen*eVecs(:,2)];
%     zAxis = [origin origin+axisLen*eVecs(:,3)];
%     plot(xAxis(1,:),xAxis(2,:),'r-','linewidth',3);
%     plot(yAxis(1,:),yAxis(2,:),'c-','linewidth',3);
%     plot(zAxis(1,:),zAxis(2,:),'y-','linewidth',3);
% end
% 
% % Draw inlier spots.
% if ~isempty(allCoordPix)
%     plot(allCoordPix(:,1),allCoordPix(:,2),'wx');
% end
% set(gca,'fontsize',20);
% hold off;
% pause(0.1);
% 
% switch job.options.debug.showPlaneFit
%     case -1
%         pause;
%     case -2
%         keyboard;
% end
% 
% end % function showPlaneFit
