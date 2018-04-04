function kitFilterCells(jobset)
% KITFILTERCELLS(EXPTS) Allows the user to remove cells whose intensities
% do not represent the general distribution from the full dataset.
%
% Copyright (c) 2018 C. A. Smith

% suppress docking
dockStatus = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal');
% suppress warnings
w = warning;
warning('off','all');

% get basic information
detC = jobset.options.coordSystemChannel;
meaC = find(jobset.options.intensity.execute);
nChans = length(meaC);
chanOrder = [2 1 3];

handles = createControls(meaC);

% get the data
job = kitLoadAllJobs(jobset);
handles.nImgs = length(job);

% get intensities over whole experiment
[ints,idxs] = getIntensities(job,meaC);
nSpots = length(ints{1});

% start progress bar
prog = kitProgress(0);

% predefine some handles required during looping
handles.imgID = 1;
handles.prevEnable = 'off';
handles.nextEnable = 'on';

% start while loop until aborted
handles.stop = 0;
while ~handles.stop
    
    % get this image ID
    iImg = handles.imgID;
    
    % check whether there is any data contained within this movie
    if ~isfield(job{iImg},'dataStruct') || ~isfield(job{iImg}.dataStruct{detC},'failed') || job{iImg}.dataStruct{detC}.failed
        handles.imgID = handles.imgID+1;
        if handles.imgID > handles.nImgs
            handles.stop = 1;
        end
        continue
    end
    
    % Change title string to represent this image number.
    handles.title.String = sprintf('1. Check image %i (of %i)',iImg,handles.nImgs);
    
    % Start waitbar.
    waitmsg = sprintf('Opening movie %i',iImg);
    hwait = waitbar(0,waitmsg);
    
    % Get current keep information.
    if isfield(job{iImg},'keep')
      handles.keep = job{iImg}.keep;
    else
      handles.keep = 1;
    end
    keeprmvCB();
    
    % Update next and prev buttons accordingly.
    handles.prevBtn.Enable = handles.prevEnable;
    handles.nextBtn.Enable = handles.nextEnable;
    
    % Show z-projected image in image panel.
    % get the image in each channel
    [~,reader] = kitOpenMovie(fullfile(job{iImg}.movieDirectory,job{iImg}.ROI.movie),'ROI');
    imgSize = kitComputeStackSize(job{iImg}.ROI.crop,job{iImg}.metadata.frameSize);
    rgbImg = zeros([imgSize(1:2) 3]);
    for iChan = 1:nChans
      img = kitReadImageStack(reader,job{iImg}.metadata,1,meaC(iChan),job{iImg}.ROI.crop,0);
      % maximum projection and contrasting
      img = max(img,[],3);
      irange(iChan,:) = stretchlim(img,[0.1 0.9999]);
      rgbImg(:,:,chanOrder(iChan)) = imadjust(img, irange(iChan,:), []);
      legend{iChan} = num2str(meaC(iChan));
    end
    set(handles.imgPanel,'cdata',rgbImg);
    clear rgbImg
    
    % Show intensity measurements in intensity tabs.
    % get the axis handle, and plot the intensities
    axes(handles.intFig);
    cla; hold on
    for iChan = 1:nChans
      % plot full dataset
      kInts = ints{iChan}(idxs~=iImg);
      xData = rand([length(kInts) 1])-0.5;
      xData = xData/2 + iChan;
      scatter(xData,kInts,10,'o','MarkerEdgeColor',[0.8 0.8 0.8]);
      % plot only this experiment
      rInts = ints{iChan}(idxs==iImg);
      xData = rand([length(rInts) 1])-0.5;
      xData = xData/2 + iChan;
      scatter(xData,rInts,10,'o','MarkerEdgeColor',[1 0 0]);
    end
    set(handles.intFig,'FontSize',8,'XTick',1:nChans);%,'YTick',[]);
    ylabel('raw intensity'); xlabel('channel');
    
    % Complete waitbar
    waitmsg = 'Finished opening movie';
    waitbar(1,hwait,waitmsg);
    close(hwait);
    
    % GUI now set up, wait for user
    uiwait(gcf);

    % back up results and save
    job{iImg}.keep = handles.keep;
    job{iImg} = kitSaveJob(job{iImg});

    % update progress
    prog = kitProgress(iImg/handles.nImgs,prog);

end

% Finish and close the GUI.
kitLog('Manual cell filtering complete');
close(gcf);

% reset docking and warning status
set(0,'DefaultFigureWindowStyle',dockStatus);
warning(w);

% Create all main controls.
function hs = createControls(meaC)
  
  % Get number of channels for intensity tabs.
  nChans = length(meaC);
    
  % Create figure.
  figw = 85;
  figh = 25;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = 'Manual cell filtering';
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';
  
  % Get pixel to character ratio.
  set(hs.fig,'Units','pixels');
  p(1,:) = get(hs.fig,'Position');
  set(hs.fig,'Units','characters');
  p(2,:) = get(hs.fig,'Position');
  hs.pixPerChar = p(1,3:4)./p(2,3:4);
  
  % Define font sizes.
  medfont = 14;
  smallfont = 12;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 1; %small horiz shift
  toplabely = figh; %top-most point
  colw = figw/2 - dx; %column width
  
  % Image check sub-title.
  labw = colw; x = dx;
  y = toplabely-lh;
  hs.title = label(hs.fig,'1. Check image',[x y labw h],medfont);
  hs.title.FontWeight = 'bold';
  
  % Create image panel.
  imgw = 42; imgh = imgw/2;
  y = 1; x = dx;
  hs.imgPanel = uicontrol(hs.fig,'Units','characters','Position',[x y imgw imgh]);

  % Intensity check sub-title.
  labw = colw; x = imgw + 2*dx;
  y = toplabely-lh;
  t = label(hs.fig,'2. Check intensities',[x y labw h],medfont);
  t.FontWeight = 'bold';
  
  % Create intensity tabs.
  panh = 15; panw = figw - (imgw+3*dx);
  y = toplabely - (panh+lh);
  hs.intPanel = uipanel(hs.fig,'Units','characters',...
      'Position',[x y panw panh],'BorderType','none');
  hs.intFig = subplot(1,1,1,'Parent',hs.intPanel);

  % Keep/ignore radio buttons.
  radw = 9; y = y-lh;
  x = imgw + 3*dx;
  hs.keepRad = uicontrol('Parent',hs.fig,'Units','characters',...
      'Style','radio','String','Keep','Position',[x y radw h],...
      'FontSize',smallfont,'Callback',@keeprmvCB);
  x = x + (radw+ddx);
  hs.rmvRad = uicontrol('Parent',hs.fig,'Units','characters',...
      'Style','radio','String','Ignore','Position',[x y radw h],...
      'FontSize',smallfont,'Callback',@keeprmvCB);
    
  % Finish, next and previous buttons
  btnw = [12 7]; btnh = 2;
  x = figw-(btnw(1)+dx); y = 1;
  hs.finishBtn = button(gcf,'Finish',[x y btnw(1) btnh],@finishCB);
  x = x-(btnw(2)+ddx);
  hs.nextBtn = button(gcf,'Next',[x y btnw(2) btnh],@nextCB);
  hs.nextBtn.Enable = 'on';
  x = x-(btnw(2)+ddx/2);
  hs.prevBtn = button(gcf,'Prev',[x y btnw(2) btnh],@prevCB);
  hs.prevBtn.Enable = 'off';
  
  movegui(hs.fig,'center');
  
end

%% Callback functions

function keeprmvCB(hObj,event)
  % update the handles
  if exist('hObj','var')
    handles.keep = strcmp(hObj.String,'Keep');
  end
  handles.keepRad.Value = handles.keep;
  handles.rmvRad.Value = ~handles.keep;
end
    
function prevCB(hObj,event)
  % update the handles
  handles.imgID = handles.imgID-1;
  if handles.imgID == 1
    handles.prevEnable = 'off';
  end
  handles.nextEnable = 'on';
  % continue the function
  uiresume(gcf);
end

function nextCB(hObj,event)
  % update the handles
  handles.imgID = handles.imgID+1;
  handles.prevEnable = 'on';
  if handles.imgID == handles.nImgs
    handles.nextEnable = 'off';
  end
  % continue the function
  uiresume(gcf);
end

function finishCB(hObj,event)
  % force stop
  handles.stop = 1;
  % continue the function
  uiresume(gcf);
end

%% Other functions

function [allInts, allIdxs] = getIntensities(jobs,chans)
    
  % get basic information
  nJobs = length(jobs);
  nChans = length(chans);
  
  % set up empty vectors
  allInts = repmat({[]},1,nChans);
  allIdxs = [];
  maxVal = nan(1,4);
  
  % loop over jobs, then channels
  for iJob = 1:nJobs
    
    % check that this job was successfully processed
    if jobs{iJob}.dataStruct{jobs{iJob}.options.coordSystemChannel}.failed
      continue
    end
      
    for iChan = 1:nChans
      
      c = chans(iChan);
      % get intensity structures
      sI = jobs{iJob}.dataStruct{c}.spotInt;
      cI = jobs{iJob}.dataStruct{c}.cellInt;
        
      % correct intensities for background  
      sI.intensity(:,iChan) = sI.intensity(:,iChan) - cI.back;
      
      % provide intensities, and compile
      ints = sI.intensity(:,iChan);
      allInts{iChan} = [allInts{iChan}; ints];
      
      % get number of spots
      nSpots = length(ints);
      
    end
    
    % provide job indices, and compile
    idxs = repmat(iJob,nSpots,1);
    allIdxs = [allIdxs; idxs];
    
  end
  
  % force intensities across channels to be over range [-Inf 1]
  for iChan = 1:nChans
    maxVal = nanmax(allInts{iChan});
    allInts{iChan} = allInts{iChan}/maxVal;
  end
  
end


end
