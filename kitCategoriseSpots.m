function kitCategoriseSpots(jobset,varargin)
% KITCATEGORISESPOTS(JOBSET,...) Displays a GUI to allow selection of spots
% from a JOBSET into a category.
%
%    KITCATEGORISESPOTS(JOBSET,...)
%
%    Options, defaults in {}:-
%
%    channel: {coordinate system channel}, 1, 2 or 3. Which channel to
%           show.
%
%    category: {'newCategory'} or string. The name under which the
%           categorisation will be saved.
%
% Copyright (c) 2018 C. A. Smith

% Process options.
opts.channel = jobset.options.coordSystemChannel;
opts.category = 'newCategory';
opts = processOptions(opts,varargin{:});

% Define colours for rectangles.
handles.col = [1 0   0;...
               0 0.75 0];
           
% Get the data.
job = kitLoadAllJobs(jobset);
handles.nMovs = length(job);
handles.chans = find(cellfun(@(x) ~strcmp(x,'none'),jobset.options.spotMode));

% Get variables from the options.
category = opts.category;
channel = opts.channel;

nFrames = job{1}.metadata.nFrames;
if nFrames > 1
  error('This function cannot yet be used for movies.')
end

% Predefine some handles required during looping
handles.movID = 1;
handles.prevEnable = 'off';
handles.nextEnable = 'on';

% Start progress.
prog = kitProgress(0);

% Start while loop until aborted
handles.stop = 0;
while ~handles.stop
        
    % get this image ID
    iMov = handles.movID;
    
    % check whether there is any data contained within this movie
    if ~isfield(job{iMov},'dataStruct') || ~isfield(job{iMov}.dataStruct{channel},'failed') || job{iMov}.dataStruct{channel}.failed
        continue
    end
    % get dataStruct
    dS = job{iMov}.dataStruct{channel};
    % get number of spots
    nSpots = length(dS.initCoord.allCoord);
    handles.nSpots = nSpots;
    
    % get categories
    if isfield(job{iMov},'categories')
        cG = job{iMov}.categories;
    else
        cG = struct;
    end
    % define keep structure
    handles.keep = zeros(nSpots,1);
    if isfield(cG,category)
        keepList = cG.(category);
        handles.keep(keepList) = 1;
    end
    
    % give nROIs to the job
    job{iMov}.nROIs = handles.nMovs;

    % show all spots - defined using tracks
    rectDims = griddedSpots(job{iMov},'channel',channel,'rawData',0);
    
    % label figure title
    figtit = sprintf('Categorise spots: Image %i%s, channel %i',job{iMov}.index,handles.nMovs,channel);
    set(gcf,'Name',figtit);
    
    % get image information
    rectPos = rectDims(:,1:2);
    rectWid = rectDims(1,3);
    handles.rectPos = rectPos;
    handles.rectWid = rectWid;
    
    % draw rectangles
    hold on
    for iSpot = 1:nSpots
        % get the colour for this spot
        keep = handles.keep(iSpot);
        icol = handles.col(keep+1,:);
        
        % draw the rectangle
        rectangle('Position',[rectPos(iSpot,:)-0.5 rectWid rectWid],...
            'EdgeColor',icol,'LineWidth',3);
    end
    
    % Buttons and labels.
    btnw = [12 7]; btnh = 2; h = 1.5;
    figpos = get(gcf,'Position');
    dx = 2.5; ddx = 1;
    % make label at top left for instructions
    x = dx; y = figpos(4)-(btnh+ddx);
    labw = 60;
    handles.instructions = label(gcf,'Click on spots to include (green) or omit (red) from the category.',[x y labw h],12);
    % add all buttons: deselect all
    x = figpos(3)-(btnw(1)+dx);
    handles.invertBtn = button(gcf,'Invert all',[x y btnw(1) btnh],@invertCB);
    x = x-(btnw(1)+ddx);
    handles.deselectBtn = button(gcf,'Deselect all',[x y btnw(1) btnh],@deselectAllCB);
    % finish, next and previous
    x = figpos(3)-(btnw(1)+dx); y = dx;
    handles.finishBtn = button(gcf,'Finish',[x y btnw(1) btnh],@finishCB);
    x = x-(btnw(2)+ddx);
    handles.nextBtn = button(gcf,'Next',[x y btnw(2) btnh],@nextMovCB);
    handles.nextBtn.Enable = handles.nextEnable;
    x = x-(btnw(2)+ddx/2);
    handles.prevBtn = button(gcf,'Prev',[x y btnw(2) btnh],@prevMovCB);
    handles.prevBtn.Enable = handles.prevEnable;
    
    % set up remove environment
    set(get(gca,'Children'),'ButtonDownFcn',@rmvCB);
    
    % GUI now set up, wait for user
    uiwait(gcf);
    
    % get the spots requiring removal
    keepList = find(handles.keep);
    
    % check for any sisters not in the list, provide warning and
    % remove if so
    incorrect = setdiff(keepList,1:nSpots);
    if ~isempty(incorrect)
        warning('The following selected spots do not exist: %s. Will ignore.',num2str(incorrect));
    end

    % process the list
    cG.(category) = keepList;
    % back up results
    job{iMov}.categories = cG;
    
    % save results
    job{iMov} = kitSaveJob(job{iMov});

    % update progress
    prog = kitProgress(iMov/handles.nMovs,prog);
    
    % close the figure for the next movie
    close(gcf);

end

kitLog('Categorisation complete.');

%% Callback functions

function rmvCB(hObj,event)
  
  % get the position of the click
  pos=get(gca,'CurrentPoint');
  xpos = pos(1,1); ypos = pos(1,2);

  % get all positions
  allPos = handles.rectPos;
  
  % get candidates using click's row position
  diffs = xpos-allPos(:,1);
  diffs(diffs<0) = NaN;
  xidx = find(diffs == nanmin(diffs));
  % get candidates using click's column position
  diffs = ypos-allPos(:,2);
  diffs(diffs<0) = NaN;
  yidx = find(diffs == nanmin(diffs));
  % get the common candidate
  idx = intersect(xidx,yidx);

  % if a click is made elsewhere, remind user how to select images
  if isempty(idx)
    handles.instructions.String = 'Click on the images to select/deselect.';
    return
  end
  
  % get the colour for this spot
  keepStat = ~handles.keep(idx);
  icol = handles.col(keepStat+1,:);

  % draw the rectangle
  rectangle('Position',[handles.rectPos(idx,:)-0.5 handles.rectWid handles.rectWid],...
      'EdgeColor',icol,'LineWidth',3);

  handles.keep(idx) = keepStat;
  
end

function invertCB(hObj,event)
  hs = handles;
  % force all stops to be ignored
  handles.keep = ~handles.keep;
  for i = 1:hs.nSpots
    
    % get the colour for this spot
    keepStat = handles.keep(i);
    jcol = handles.col(keepStat+1,:);
    
    % draw the rectangle
    rectangle('Position',[hs.rectPos(i,:)-0.5 hs.rectWid hs.rectWid],...
        'EdgeColor',jcol,'LineWidth',3);
  end
end

function deselectAllCB(hObj, event)
  hs = handles;
  % force all stops to be ignored
  handles.keep = 0.*handles.keep;
  for i = 1:hs.nSpots

    % get the colour for this spot
    keepStat = handles.keep(i);
    jcol = handles.col(keepStat+1,:);

    % draw the rectangle
    rectangle('Position',[hs.rectPos(i,:)-0.5 hs.rectWid hs.rectWid],...
        'EdgeColor',jcol,'LineWidth',3);
  end
end

function prevMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID-1;
  if handles.movID == 1
    handles.prevEnable = 'off';
  end
  handles.nextEnable = 'on';
  handles.nextChan = 0;
  % continue the function
  uiresume(gcf);
end

function nextMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID+1;
  handles.prevEnable = 'on';
  if handles.movID == handles.nMovs
    handles.nextEnable = 'off';
  end
  handles.nextChan = 0;
  % continue the function
  uiresume(gcf);
end

function finishCB(hObj,event)
  % force stop
  handles.stop = 1;
  % continue the function
  uiresume(gcf);
end

end
