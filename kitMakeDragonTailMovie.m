function kitMakeDragonTailMovie(job,varargin)
% KITMAKEDRAGONTAILMOVIE(JOB,...) Makes a .gif file showing dragon tails
% for all tracks in a given JOB.
%
% Copyright (c) 2017 C. A. Smith

% default options
opts.contrast = {[0.1 1]};
opts.delay = 1; % in s
opts.filename = 'dragonTailMovie.gif';
opts.scaleBarLabel = 0;
opts.scaleBarSize = 3; % in µm, length of scale bar
opts.timePoints = [];
opts.tailLength = 8;
opts.timeStamp = 1;
opts.topRight = [];
opts.bottomRight = [];
% process user-defined options
opts = processOptions(opts,varargin{:});

% get image size information
if isempty(opts.timePoints)
  opts.timePoints = 1:job.metadata.nFrames;
end
imageSize = job.ROI(job.index).cropSize(1:2);
timeLapse = job.metadata.frameTime(1,2); %in s
if opts.timeStamp
  % check size of movie to deduce style of timeStamp
  hourStamp = (timeLapse*max(opts.timePoints))>3600;
end

% get length of scale bar
pixelSize = job.metadata.pixelSize(1);
barLength = opts.scaleBarSize/pixelSize(1);

for iFrame = opts.timePoints
  
  % plot the image with dragon tails
  kitShowDragonTails(job,'timePoint',iFrame,'tailLength',opts.tailLength,'contrast',opts.contrast);
  
  % overlay requested information
  indent = 10;
  if barLength > 0
    line([indent indent+barLength], [imageSize(1)-indent imageSize(1)-indent],...
        'Color','w','LineWidth',5);
    if opts.scaleBarLabel
      scaleLabel = [num2str(opts.scaleBarSize) ' µm'];
      text(indent+(barLength/2),imageSize(1)-indent-5,scaleLabel,...
          'Color','w','HorizontalAlignment','center','FontSize',20);
    end
  end
  if opts.timeStamp
    % get time in days, and convert
    secs = (iFrame-1)*timeLapse;
    secs = secs/(24*60*60);
    timeLabel = datestr(secs,'HH:MM:SS');
    if ~hourStamp
      timeLabel = timeLabel(4:end);
    end
    text(indent,indent,timeLabel,'Color','w','FontSize',30);
  end
  if ~isempty(opts.topRight)
    % print label requested in top right
    if iscell(opts.topRight)
      vertIndent = indent+((length(opts.topRight)-1)*5);
    else
      vertIndent = indent;
    end
    text(imageSize(2)-indent,vertIndent,opts.topRight,...
        'Color','w','FontSize',25,'HorizontalAlignment','right');
  end
  if ~isempty(opts.bottomRight)
    % print label requested in top right
    text(imageSize(2)-indent,imageSize(1)-indent,opts.bottomRight,...
        'Color','w','FontSize',25,'HorizontalAlignment','right');
  end
  
  % draw the figure as an image
  drawnow
  frame = getframe(1);
  outputImg{iFrame} = frame2im(frame);
  % convert it to indices
  [A,map] = rgb2ind(outputImg{iFrame},256);
  
  if iFrame == opts.timePoints(1)
    % write first timepoint
    imwrite(A,map,opts.filename,'gif','LoopCount',Inf,'DelayTime',opts.delay);
    
  else
    % append later timepoints to the current movie
    imwrite(A,map,opts.filename,'gif','WriteMode','append','DelayTime',opts.delay);
    
  end
  
end





end