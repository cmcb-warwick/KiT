function kitVisualiseTrackedSubset(job, channel, subset, makeMovie)
%
% plot a subset of trajectories on the original movie
%example: kitPlotTracks(job,'subset',[1,2,4,6,10,216,246],'faintSubset',1:332)
%Jonathan U. Harrison 2019-06-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    makeMovie=0;
end

if nargin<3
    subset=[];
end

if nargin<2
    channel=1;
end

if isfield(job,'metadata')
    if iscell(job.metadata)
        job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
end


nTimePoints = job.metadata.nFrames;

coords = cell(nTimePoints,1);
for iTrack=subset
    for iFrame=1:nTimePoints
        %TODO: fix use of temporary object by fixing coord conversion size
        %issue
        try
          coordsTMP = kitCoordsToImageCoords(job,channel,...
            job.dataStruct{channel}.trackList(iTrack).coords(iFrame,1:3),iFrame);
        catch
          coordsTMP = NaN(1,3);
        end
    coords{iFrame}(iTrack,1:3) =  coordsTMP(1,:);    
    end
end
  if makeMovie
    h = figure;
    if ~verLessThan('matlab','9.5') %this property only added for matlab R2018a onwards
      h.WindowState = 'fullscreen'; %make figure full screen 
    end
    outname = kitGenerateOutputFilename(job);
    vidfile = VideoWriter(sprintf('%sTrackedSubsetMovie.mp4',outname),'MPEG-4');
    open(vidfile);
  end
  for iFrame = 1:nTimePoints
    showSpots(job,reader,channel,iFrame,coords{iFrame})
    if makeMovie
      drawnow;
      F(iFrame) = getframe(h);
      writeVideo(vidfile,F(iFrame));
    end
  end
  if makeMovie
    close(vidfile)
  end

end

function showSpots(job,reader,channel,frameNum,allCoordPix)
% Based on debugging visualization from kitFitPlane

img = kitReadImageStack(reader, job.metadata, frameNum, channel, job.ROI.crop);
% Max project image.
img = max(img,[],3);

% Draw image.
imshow(img,[]);
title(sprintf('Plane fit frame %d/%d (x red, y cyan, z yellow)',...
              frameNum,job.metadata.nFrames));
hold on;

% Draw spots.
if ~isempty(allCoordPix)
  plot(allCoordPix(:,1),allCoordPix(:,2),'wx');
end
set(gca,'fontsize',20);
hold off;
pause(0.1);

end 
