function kitMultiJobsetRun(jobset)
% KITMULTIJOBSETRUN Display GUI to enable running spot detection of
% multiple jobs.
%
% Copyright (c) 2017 C. A. Smith

% Check whether user has asked for help.
if nargin<1
  jobset=[];
end

% Setup GUI.
dockStatus = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal');
handles = createControls(jobset);
handles.jobsets{1} = jobset;
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);
set(0,'DefaultFigureWindowStyle',dockStatus);

%% NESTED FUNCTIONS
function hs = createControls(jobset)
  
  % Create figure.
  figw = 50;
  figh = 20;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = 'Run multiple jobsets';
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  % Define font sizes.
  medfont = 14;
  smallfont = 12;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 0.5; %small horizontal shift
  toplabely = figh; %top-most point
  
  % Set up initial positions.
  x = dx;
  w = figw-2*dx;
  y = toplabely-lh;
  
  % Running tools sub-title.
  labw = w;
  t = label(hs.fig,'1. List multiple jobsets',[x y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Buttons for adding and removing runs.
  btnw = 12; btnh = 2;
  hs.newBtn = button(hs.fig,'New run',[x y btnw btnh],@newRunCB);
  x = x+(btnw+ddx);
  hs.loadBtn = button(hs.fig,'Load run',[x y btnw btnh],@loadRunCB);
  x = x+(btnw+ddx);
  hs.rmvBtn = button(hs.fig,'Remove',[x y btnw btnh],@rmvRunCB);
  
  % Listing loaded jobsets.
  panw = w; panh = 10;
  y = y-panh; x = dx;
  hs.jobsetList = uicontrol(hs.fig,'Style','listbox','Units','characters',...
      'Position',[x y panw panh],'Max',Inf,'Min',0,'FontSize',smallfont);
  jScrollPane = findjobj(hs.jobsetList); % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane.setHorizontalScrollBarPolicy(32);
  if ~isempty(jobset)
    hs.jobsets{1} = jobset;
    filename = jobset.filename;
    maxMovLen = 55;
    filename = strshorten(filename,maxMovLen);
    hs.jobsetList.String{1} = filename;
  end
  
  % Run button.
  x = dx; y = 1+lh;
  % Running tools sub-title.
  labw = w;
  t = label(hs.fig,'2. Run all jobs in the list',[x y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-lh;
  hs.runBtn = button(hs.fig,'Run jobs',[x y btnw btnh],@executeCB,smallfont);
  
  % Close button.
  x = figw-(btnw+dx); y = 1;
  hs.closeBtn = button(hs.fig,'Close',[x y btnw btnh],@closeCB);

  movegui(hs.fig,'center');
  
end

function newRunCB(hObj,event)
  kitLog('Opening jobset setup window')
  hs = handles;
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  pause(1);
  
  % Small GUI to ask user which type of job to set up.
  mode = chooseMode;
  
  jobset = kitSetupJob(mode);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
  % Check whether user cancelled setup.
  if isfield(jobset,'cancel')
    return
  end
  
  % push jobset name to label
  filename = jobset.filename;
  maxMovLen = 55;
  filename = strshorten(filename,maxMovLen);
  hs.jobsetList.String{end+1} = filename;
  % back up all information
  hs.jobsets{end+1} = jobset;
  handles = hs;
end

function loadRunCB(hObj,event)
  kitLog('Loading jobset')
  hs = handles;
  % run GUI to make new jobset
  [filename,pathname] = uigetfile('*.mat','Select a Sid jobset file');
  jobset = kitLoadJobset(fullfile(pathname,filename));
  % push jobset name to label
  filename = jobset.filename;
  maxMovLen = 55;
  filename = strshorten(filename,maxMovLen);
  hs.jobsetList.String{end+1} = filename;
  hs.runBtn.String = 'Run';
  % back up all information
  if ~isfield(hs,'jobsets')
    hs.jobsets{1} = jobset;
  else
    hs.jobsets{end+1} = jobset;
  end
  handles = hs;
end

function rmvRunCB(hObj,event)
  hs = handles;
  % find which jobsets are highlighted
  rmvidx = hs.jobsetList.Value;
  if isempty(rmvidx)
    errorbox('Must first select jobset(s) to remove');
    return
  end
  jobsetList = hs.jobsetList.String;
  if length(rmvidx) == length(hs.jobsets)
    jobsetList = '';
    % remove selected jobsets
    hs.jobsets = [];
  else
    jobsetList(rmvidx) = [];
    set(hs.jobsetList,'Value',1);
    % remove selected jobsets
    hs.jobsets(rmvidx) = [];
  end
  set(hs.jobsetList,'String',jobsetList);
  handles = hs;
end

function executeCB(hObj,event)
  if ~isfield(handles,'jobsets') || isempty(handles.jobsets)
    errorbox('No jobs yet created or loaded.')
    return
  end
  
  % Turn off the GUI during processing.
  handles.runBtn.String = 'Running...';
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  pause(1);
  
  % Run the job.
  kitRunAllJobs(handles.jobsets);
  
  % We turn back on the interface
  set(guiObj,'Enable','on');
  handles.runBtn.String = 'Re-run';
end

function closeCB(hObj,event)
  uiresume(gcf);
end

%% Other functions

function mode = chooseMode()
  
  % List all modes.
  allModes = {'Tracking','Single timepoint','Chromatic shift'};
  allModesJS = {'zandt','zonly','chrshift'};
    
  % Create figure.
  figw = 50;
  figh = 5;
  f = figure('Resize','off','Units','characters','Position',[100 35 figw figh]);
  f.DockControls = 'off';
  f.MenuBar = 'none';
  f.Name = 'Choose a job type';
  f.NumberTitle = 'off';
  f.IntegerHandle = 'off';
  f.ToolBar = 'none';
  
  % Define font sizes.
  smallfont = 12;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 0.5; %small horizontal shift
  toplabely = figh; %top-most point
  
  % Set up initial positions.
  x = dx;
  w = figw-2*dx;
  y = toplabely-lh;
  
  % Print choices.
  labw = w;
  label(f,'Choose which type of job you wish to set up:',[x y labw h],smallfont);
  y = y-lh;
  btnw = (w-2*ddx)/3; btnh = 2;
  for i=1:3
    button(f,allModes{i},[x y btnw btnh],@chooseModeCB,smallfont);
    x = x+(btnw+ddx);
  end
  
  movegui(f,'center');
  uiwait(f);
  
  % Get the mode for the jobset.
  mode = allModesJS{mode};
    
  function chooseModeCB(hObj,event)
    mode = find(cellfun(@(x) strcmp(x,hObj.String),allModes));
    uiresume(f);
    close(f);
  end
  
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

end % kitSetupJob