function spotSel = kitSelectData(jobset,varargin)
% KITSELECTDATA(JOBSET,...) Displays a GUI to allow removal of erroneous
% spots from a JOBSET.
%
% SPOTS MUST BE PAIRED PRIOR TO THIS FILTERING!!!
%
% Copyright (c) 2018 C. A. Smith

% Define colours for rectangles.
handles.col = [1 0   0;...
               0 0.75 0];
           
% Get the data.
job = kitLoadAllJobs(jobset);
handles.nMovs = length(job);
handles.chans = find(cellfun(@(x) ~strcmp(x,'none'),jobset.options.spotMode));

% Get basic information from the jobset.
opts = jobset.options;
handles.chanID = opts.coordSystemChannel;

nFrames = job{1}.metadata.nFrames;
if nFrames > 1
  error('This function cannot yet be used for movies.')
end

% Predefine some handles required during looping
handles.movID = 1;
handles.prevEnable = 'off';
handles.nextEnable = 'on';
allSels = [];
trackSels = [];

% Start progress.
prog = kitProgress(0);

% Start while loop until aborted
handles.stop = 0;
while ~handles.stop
        
    if handles.movID > handles.nMovs
      handles.stop=1;
      continue
    end
    
    % get this image ID
    iMov = handles.movID;
    % get this channel ID
    iChan = handles.chanID;

    % check whether there is any data contained within this movie
    if ~isfield(job{iMov},'dataStruct') || length(job{iMov}.dataStruct)<iChan || ...
            ~isfield(job{iMov}.dataStruct{iChan},'failed') || job{iMov}.dataStruct{iChan}.failed
        handles.movID = handles.movID+1;
        if handles.movID > handles.nMovs
          handles.stop = 1;
        end
        continue
    end
    % get dataStruct
    dS = job{iMov}.dataStruct{iChan};
    % get the full initCoord and spotInt
    iC = dS.initCoord;
    % get IDs for all spots not previously filtered
    nonNaNs = find(~isnan(iC.allCoord(:,1)));
    
    % back up the full initCoord and spotInt
    if isfield(dS,'rawData')
      iC = dS.rawData.initCoord;
      if isfield(dS.rawData,'spotInt')
        sI = dS.rawData.spotInt;
      end
    elseif isfield(dS,'spotInt')
        sI = dS.spotInt;
    end
    raw.initCoord = iC;
    if exist('sI','var')
        raw.spotInt = sI;
    end
    dS.rawData = raw;
    job{iMov}.dataStruct{iChan} = dS;
    
    % get number of spots
    nSpots = size(dS.initCoord.allCoord,1);
    handles.nSpots = nSpots;

    % show all spots - defined using tracks
    rectDims = griddedSpots(job{iMov},'channel',handles.chanID);
    
    % get image information
    rectPos = rectDims(:,1:2);
    rectWid = rectDims(1,3);
    handles.rectPos = rectPos;
    handles.rectWid = rectWid;
    handles.keep = ismember(1:nSpots,nonNaNs);
    
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
    handles.instructions = label(gcf,'Click on spots to keep (green) or remove (red).',[x y labw h],12);
    % add all buttons: finish, next and previous
    x = figpos(3)-(btnw(1)+dx);
    handles.nextBtn = button(gcf,'Next',[x y btnw(2) btnh],@nextMovCB);
    handles.nextBtn.Enable = handles.nextEnable;
    % invert all
    x = figpos(3)-(btnw(1)+dx); y = dx;
    handles.invertBtn = button(gcf,'Invert all',[x y btnw(1) btnh],@invertCB);
    
    % set up remove environment
    set(get(gca,'Children'),'ButtonDownFcn',@rmvCB);
    
    % GUI now set up, wait for user
    uiwait(gcf);
    
    % get the spots requiring removal
    keepList = find(handles.keep);

    % compile selection information
    thisSel = [];
    thisSel(:,1) = repmat(iMov,length(keepList),1);
    thisSel(:,2) = keepList(:);
    allSels = [allSels; thisSel];
    allSels = unique(allSels,'rows');
    
    % convert to trackIDs
    allSpotIDs = cat(1,job{iMov}.dataStruct{iChan}.trackList.featIndx);
    keepList_trk = find(ismember(allSpotIDs,keepList(:)));
    trackSels = [trackSels; [repmat(iMov,length(keepList_trk),1) keepList_trk(:)] ];
    
    % update progress
    prog = kitProgress(iMov/handles.nMovs,prog);
    
    % close the figure for the next movie
    close(gcf);

end

% Store results in structure.
spotSel.dataType = 'initCoord';
spotSel.rawSelection{1} = allSels;
spotSel.selection{1} = trackSels;

kitLog('Manual filtering complete.');

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

function nextMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID+1;
  if handles.movID == handles.nMovs
      handles.nextBtn.String = 'Finish';
  end
  handles.nextChan = 0;
  % continue the function
  uiresume(gcf);
end

end


    