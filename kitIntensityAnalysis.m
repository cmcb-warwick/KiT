function varargout=kitIntensityAnalysis()
% KITINTENSITYANALYSIS Displays GUI to allow user to run analysis of
% pre-processed intensity measurements. Either jobsets or intrameasurement
% structures can be used.
%
% Copyright (c) 2017 C. A. Smith

% Set up GUI.
dockStatus = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal');
handles = createControls();
updateControls();
handles.fig.Visible = 'on';
uiwait(gcf);
if nargout==1
    varargout{1} = analysis;
end
set(0,'DefaultFigureWindowStyle',dockStatus);

%% NESTED FUNCTIONS

% Create all main controls.
function hs = createControls()
  
  % Create figure.
  figw = 90;
  figh = 34;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = 'Intensity analysis';
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
  tinyfont = 10;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 1; %small horiz shift
  toplabely = figh; %top-most point
  colw = figw/2 - dx; %column width
  
  % Conditions sub-title.
  labw = colw; x = dx;
  y = toplabely-lh;
  t = label(hs.fig,'1. Select data per condition',[x y labw h],medfont);
  t.FontWeight = 'bold';
  
  % Data type selection.
  txtw = 10; radw = 17.5;
  x = dx;
  y = y-lh;
  hs.condsText = label(hs.fig,'Data type:',[x y txtw h],smallfont);
  x = x+(txtw+ddx);
  hs.dataType{1} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String','intra-measurements','Position',[x y radw h],'Callback',@dataTypeCB);
  x = x+(radw+ddx);
  hs.dataType{2} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String','jobsets','Position',[x y radw h],'Callback',@dataTypeCB,'Value',1);
  
  % Add/delete conditions.
  btnw = 6; txtw = 12.5;
  x = dx;
  y = y-lh;
  hs.condsText = label(hs.fig,'Conditions:',[x y txtw h],smallfont);
  hs.condsText.FontWeight = 'bold';
  x = x+(txtw);
  hs.addcond = button(hs.fig,'Add',[x y btnw h],@addcondCB,tinyfont);
  x = x+(btnw+dx/4);
  hs.delcond = button(hs.fig,'Delete',[x y btnw h],@delcondCB,tinyfont);
  
  % Create condition tabs.
  panw = (figw-dx)/2; panh = 20;
  pany = y-(panh+dx/4);
  x = dx;
  hs.tabPanel = uipanel(hs.fig,'Units','characters',...
      'Position',[x pany panw panh],'FontSize',smallfont,'BorderType','none');
  hs.tabP = uitabgroup('Parent', hs.tabPanel);
  hs.tabs{1} = uitab('Parent', hs.tabP, 'Title', '1');
  hs = createTabControls(1,hs);
  
  % Options sub-title.
  colx = (figw+dx)/2;
  x = colx+dx; labw = colw;
  y = toplabely-lh;
  t = label(hs.fig,'2. Options',[x y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Control channel.
  nChans = 3;
  txtw = 15;
  label(hs.fig,'Control channel',[x y txtw h],smallfont);
  radw = 5; x = figw-(radw-dx);
  for i=1:nChans
    hs.controlCh{i} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String',num2str(i),'Position',[x-radw*(4-i) y radw h],'Callback',@contChCB);
  end
  hs.controlCh{1}.Value = 1;
  hs.controlChNum = 1;
  % Channels to analyse.
  y = y-h;
  x = colx+dx; txtw = 15;
  label(hs.fig,'To analyse...',[x y txtw h],smallfont);
  tickw = 5; x = figw-(tickw+dx);
  for i=nChans:-1:1
    hs.analyseCh{i} = checkbox(hs.fig,num2str(i),[x y tickw h],@analyseChCB,tinyfont);
    x = x-tickw;
  end
  hs.analyseCh{1}.Value = 1;
  hs.analyseCh{1}.Enable = 'off';
  % Channel labels.
  y = y-h;
  x = colx+dx; txtw = 15;
  label(hs.fig,'Channel labels',[x y txtw h],smallfont);
  x = x+txtw+2*dx; txtw = 5; editw = figw-(x+txtw+dx);
  for i=1:nChans
    hs.chlabelsTxt{i} = label(hs.fig,['ch ' num2str(i)],[x y txtw h],tinyfont);
    hs.chlabels{i} = editbox(hs.fig,'',[x+txtw y editw h],tinyfont);
    y = y-h;
  end
  
  % Plot style buttons.
  x = colx+dx; txtw = 15;
  label(hs.fig,'Select plot style',[x y txtw h],smallfont);
  y = y-h;
  btnw = ((figw-colx)-(2*ddx+2*dx))/3;
  btnh = btnw*0.45;
  y = y-(btnh-h);
  hs.boxWhiskBtn = uicontrol('style','pushbutton','Units','characters',...
      'position',[x y btnw btnh],'Callback',@boxWhiskBtnCB);
  hs.boxWhiskJava = java(findjobj(hs.boxWhiskBtn));
  changeBorder(hs.boxWhiskJava,1);
  hs.plotStyle = 'boxwhisk';
  [I,~]=imread('private/boxWhiskBtn.png','BackgroundColor',[1 1 1]);
  I=imresize(I,0.8*[btnw btnh].*hs.pixPerChar);
  set(hs.boxWhiskBtn,'cdata',I);
  x = x+btnw+ddx;
  hs.barBtn = uicontrol('style','pushbutton','Units','characters',...
      'position',[x y btnw btnh],'Callback',@barBtnCB);
  hs.barJava = java(findjobj(hs.barBtn));
  changeBorder(hs.barJava,0);
  [I,~]=imread('private/barBtn.png','BackgroundColor',[1 1 1]);
  I=imresize(I,0.8*[btnw btnh].*hs.pixPerChar);
  set(hs.barBtn,'cdata',I);
  x = x+btnw+ddx;
  hs.scatterBtn = uicontrol('style','pushbutton','Units','characters',...
      'position',[x y btnw btnh],'Callback',@scatterBtnCB);
  hs.scatterJava = java(findjobj(hs.scatterBtn));
  changeBorder(hs.scatterJava,0);
  [I,~]=imread('private/scatterBtn.png','BackgroundColor',[1 1 1]);
  I=imresize(I,0.8*[btnw btnh].*hs.pixPerChar);
  set(hs.scatterBtn,'cdata',I);
  
  % Figure option.
  y = y-lh;
  txtw = 30;
  x = figw-txtw;
  hs.closeFigs = checkbox(hs.fig,'Leave figures open afterwards',[x y txtw h],[],tinyfont);
  hs.closeFigs.Value = 0;
  
  % Running sub-title.
  labw = colw; x = dx;
  y = 1+2*lh;
  t = label(hs.fig,'3. Saving and running',[x y labw h],medfont);
  t.FontWeight = 'bold';
  
  % Save directory.
  x = dx; btnw = 15; y = 1+lh;
  button(hs.fig,'Save directory:',[x y btnw h],@savedirCB);
  x = x+btnw+ddx; txtw = 27.5;
  hs.savedirText = label(hs.fig,'./',[x y txtw h],smallfont);
  hs.savedirText.HorizontalAlignment = 'right';
  x = x+txtw;
  btnw=10; txtw = figw-(x+ddx+btnw+dx); %btnw here refers to the 'Run' button
  hs.filename = editbox(hs.fig,'filename',[x y txtw h],smallfont);
  
  % 'Run' button.
  btnw = 10; btnh = 2;
  x = figw-(btnw+dx);
  hs.runbtn = button(hs.fig,'Run',[x y btnw btnh],@runCB,smallfont);
  
  % 'Close' button.
  btnw = 10; btnh = 2;
  x = figw-(btnw+dx);
  y = y-lh;
  hs.runbtn = button(hs.fig,'Close',[x y btnw btnh],@closeCB,smallfont);
  
  movegui(hs.fig,'center');
  
end
  
% Create each condition tab.
function hs = createTabControls(tabidx,hs)
  
  % Get tab handle.
  tabhs = hs.tabs{tabidx};
  panw = hs.tabPanel.Position(3);
  panh = hs.tabPanel.Position(4);
    
  % Define font sizes.
  smallfont = 12;
  tinyfont = 10;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  toplabely = panh-2.5*h; %top-most point
  w = panw-3*dx;
  
  % Condition label.
  x = dx;
  txtw = 15;
  y = toplabely-h;
  label(tabhs,'Condition label',[x y txtw h],smallfont);
  x = x+(txtw+dx);
  txtw = 10;
  hs.condlabel{tabidx} = editbox(tabhs,'',[x y txtw h],smallfont);
  x = dx;
  tickw = panw;
  y = y-h;
  hs.normConds{tabidx} = checkbox(tabhs,'Normalise to this condition',[x y tickw h],@normcondCB,tinyfont);
  
  % Populate jobset box.
  x = dx; btnw = 15;
  y = y-lh;
  button(tabhs,'Find jobset(s)',[x y btnw h],@jobsetdirCB,smallfont);
  x = x+(btnw+dx); labw = panw-(btnw+3*dx);
  hs.allfiles{tabidx} = checkbox(tabhs,'Show all .mat files',[x y labw h],[],tinyfont);
  x = dx; y = 1;
  hs.jobsets{tabidx} = uicontrol(tabhs,'Style','listbox','Units','characters','Position',[x y w panh/2],'Max',Inf,'Min',0);
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.jobsets{tabidx});
  jScrollPane.setHorizontalScrollBarPolicy(32);
    
end

% Update control status based on contents of jobset.
function updateControls()
  
  % list all required controls here
  contChCB();
  
end

%% Check controls for consistent input.
function tf=checkControls()
  hs = handles;
  
  v = str2double(hs.autoRadiidt.String);
  if (hs.autoRadii.Value == 1) && (~isfinite(v) || v<=0)
    errorbox('Invalid value for frame dt. Should be a positive number.')
    tf = false;
    return
  end

  v = str2double(hs.autoRadiiAvgDisp.String);
  if (hs.autoRadii.Value == 1) && (~isfinite(v) || v<0)
    errorbox('Invalid value for avg. disp. of spots. Should be a positive number.')
    tf = false;
    return
  end

  v1 = str2double(hs.minSearchRadius.String);
  v2 = str2double(hs.maxSearchRadius.String);
  if (hs.autoRadii.Value ~= 1) && (~isfinite(v1) || v1 > v2 || v1 < 0)
    errorbox('Invalid value min search radius. Should be a positive number less than max search radius.')
    tf = false;
    return
  end
  if (hs.autoRadii.Value ~= 1) && (~isfinite(v2) || v2 < v1 || v2 < 0)
    errorbox('Invalid value max search radius. Should be a positive number less than min search radius.')
    tf = false;
    return
  end

  v = str2double(hs.maxSisterAlignmentAngle.String);
  if (hs.useSisterAlignment.Value == 1) && (~isfinite(v) || v < 0 || v > 180)
    errorbox('Invalid value for max angle between sisters and plate normal. Should be a positive number.')
    tf = false;
    return
  end

  v = str2double(hs.maxSisterDist.String);
  if ~isfinite(v) || v < 0
    errorbox('Invalid value for max intersister distance. Should be a positive number.')
    tf = false;
    return
  end

  v = str2double(hs.minSisterTrackOverlap.String);
  if ~isfinite(v) || v < 0
    errorbox('Invalid value for min overlap between sister tracks. Should be a positive number.')
    tf = false;
    return
  end

  v = str2double(hs.minSpotsPerFrame.String);
  if ~isfinite(v) || v < 0
    errorbox('Invalid value for min spots per frame. Should be a positive number.')
    tf = false;
    return
  end

  v = str2double(hs.adaptiveLambda.String);
  if (~isfinite(v) || v < 0)
    errorbox('Invalid value for spot weight. Should be a positive number or zero.')
    tf = false;
    return
  end

  v = str2double(hs.maxMmfTime.String);
  if (hs.mmfAddSpots.Value == 1) && (~isfinite(v) || v < 0)
    errorbox('Invalid value for min spots per frame. Should be a positive number or zero.')
    tf = false;
    return
  end

  tf = true;
end

%% All callback functions.
function dataTypeCB(hObj,event)
    
  if exist('hObj','var')  
    hs = handles;
    type = hObj.String;
    switch type
      case 'intra-measurements'
        hs.dataType{1}.Value = 1;
        hs.dataType{2}.Value = 0;
        % Enable channels 1 and 2, and change labels.
        hs.controlChNum = 1;
        hs.controlCh{1}.Value = 1;
        hs.controlCh{1}.String = 'i';
        hs.analyseCh{1}.Value = 1;
        hs.analyseCh{1}.String = 'i';
        hs.chlabels{1}.Enable = 'on';
        hs.chlabelsTxt{1}.String = 'inner';
        hs.chlabelsTxt{1}.Enable = 'on';
        hs.controlCh{2}.Value = 0;
        hs.controlCh{2}.String = 'o';
        hs.analyseCh{2}.Value = 1;
        hs.analyseCh{2}.String = 'o';
        hs.chlabels{2}.Enable = 'on';
        hs.chlabelsTxt{2}.String = 'outer';
        hs.chlabelsTxt{2}.Enable = 'on';
        % silence channel 3
        hs.controlCh{3}.Value = 0;
        hs.controlCh{3}.String = '';
        hs.controlCh{3}.Enable = 'off';
        hs.analyseCh{3}.Value = 0;
        hs.analyseCh{3}.String = '';
        hs.analyseCh{3}.Enable = 'off';
        hs.chlabels{3}.Enable = 'off';
        hs.chlabelsTxt{3}.String = '';
        hs.chlabelsTxt{3}.Enable = 'off';
      case 'jobsets'
        hs.dataType{1}.Value = 0;
        hs.dataType{2}.Value = 1;
        % Enable channels 1 and 2, and change labels.
        hs.controlChNum = 1;
        hs.controlCh{1}.Value = 1;
        hs.analyseCh{1}.Value = 1;
        hs.analyseCh{1}.Enable = 'off';
        hs.controlCh{3}.Enable = 'on';
        for notChan = 2:3
          hs.controlCh{notChan}.Value = 0;
          hs.analyseCh{notChan}.Enable = 'on';
          hs.analyseCh{notChan}.Value = 0;
        end
        for iChan = 1:3
          hs.chlabelsTxt{iChan}.String = sprintf('Ch. %i',iChan);
          hs.controlCh{iChan}.String = num2str(iChan);
          hs.analyseCh{iChan}.String = num2str(iChan);
        end
        analyseChCB();
    end
    
    
    handles = hs;
  end
    
end

function savedirCB(hObj,event)
  
  % ask the user to select a directory
  dirname = uigetdir([],'Select directory for save analysis');
  % add a '/' to the end, then shorten if needs be
  if ~isempty(dirname)
    % save the full directory for later
    handles.savedir = dirname;
    % shorten for the jobset
    dirname = [dirname filesep];
    maxDirLen = 30;
    dirname = strshorten(dirname,maxDirLen);
    % push this to the text on the left of the filename
    set(handles.savedirText,'String',dirname);
  end
  
end

function normcondCB(hObj,event)
    
    % get number of conditions and selected condition ID
    nConds = length(handles.tabP.Children);
    condtit = strrep(handles.tabP.SelectedTab.Title,'*','');
    condidx = str2double(condtit);
    % adjust ticks and aesthetics
    if handles.normConds{condidx}.Value
        handles.tabs{condidx}.Title = [handles.tabs{condidx}.Title '*'];
        for notCond = setdiff(1:nConds,condidx)
            handles.normConds{notCond}.Value = 0;
            handles.tabs{notCond}.Title = num2str(notCond);
        end
    else
        handles.tabs{condidx}.Title = num2str(condidx);
    end
    
end

function addcondCB(hObj,event)
    hs = handles;
    
    % get the number of conditions, and set new tab number
    nConds = length(hs.tabP.Children); condidx = nConds+1;
    % add the tab and switch to this one, then create the controls
    hs.tabs{condidx} = uitab('Parent', hs.tabP, 'Title', num2str(condidx));
    hs = createTabControls(condidx,hs);
    hs.jobsetdir{condidx} = [];
    
    pause(0.5); % required to avoid crashing
    hs.tabP.SelectedTab = hs.tabs{condidx};
    handles = hs;
end

function delcondCB(hObj,event)
    % get handles
    hs = handles;
    % get the number of conditions and selected tab
    nConds = length(hs.tabP.Children);
    condtit = strrep(hs.tabP.SelectedTab.Title,'*','');
    condidx = str2double(condtit);
    % delete the tab, and relabel all tabs to be 1:nConds
	delete(hs.tabs{condidx});
    hs.tabs(condidx) = [];
    hs.condlabel(condidx) = [];
    hs.normConds(condidx) = [];
    hs.jobsets(condidx) = [];
    hs.jobsetdir(condidx) = [];
    hs.allfiles(condidx) = [];
    for i=1:nConds-1
       hs.tabs{i}.Title = num2str(i);
    end
    pause(0.5); % required to avoid crashing
    % give back handles
    handles = hs;
end

function jobsetdirCB(hObj,event)
    
    % ask user to select a diretory
    dirname = uigetdir([],'Select directory containing jobsets for this condition');
    if ~isempty(dirname)
        % get the selected tab ID and fill the jobset box with .mat files
        condtit = strrep(handles.tabP.SelectedTab.Title,'*','');
        condidx = str2double(condtit);
        handles.jobsetdir{condidx}.String = dirname;
        populateJSbox(condidx);
    end
    
end

function contChCB(hObj,event)
  
  if exist('hObj','var')
      
    hs = handles;
    
    % set new control channel
    chan = str2double(hObj.String);
    hs.controlChNum = chan;
    hs.analyseCh{chan}.Value = 1;
    hs.analyseCh{chan}.Enable = 'off';
    
    % adjust control and analysis channel information accordingly
    if hs.dataType{1}.Value
      notChan = setdiff(1:2,chan);
      hs.controlCh{notChan}.Value = 0;
      hs.analyseCh{notChan}.Enable = 'on';
      hs.analyseCh{notChan}.Value = 1;
    else  
      for notChan = setdiff(1:3,chan)
        hs.controlCh{notChan}.Value = 0;
        hs.analyseCh{notChan}.Enable = 'on';
        hs.analyseCh{notChan}.Value = 0;
      end
      analyseChCB();
    end
    
    handles = hs;
    
  end
  
end

function analyseChCB(hObj,event)
    for iChan = 1:3
        if handles.analyseCh{iChan}.Value
            handles.chlabels{iChan}.Enable = 'on';
            handles.chlabelsTxt{iChan}.Enable = 'on';
        else
            handles.chlabels{iChan}.String = '';
            handles.chlabels{iChan}.Enable = 'off';
            handles.chlabelsTxt{iChan}.Enable = 'off';
        end
    end
end

function boxWhiskBtnCB(hObj,event)
  % check whether plot style is already boxWhisker
  if strcmp(handles.plotStyle,'boxwhisk')
    return
  end
  % reverse border on each
  changeBorder(handles.boxWhiskJava,1);
  changeBorder(handles.barJava,0);
  changeBorder(handles.scatterJava,0);
  % change plot style
  handles.plotStyle = 'boxwhisk';
end

function barBtnCB(hObj,event)
  % check whether plot style is already boxWhisker
  if strcmp(handles.plotStyle,'bars')
    return
  end
  % reverse border on each
  changeBorder(handles.boxWhiskJava,0);
  changeBorder(handles.barJava,1);
  changeBorder(handles.scatterJava,0);
  % change plot style
  handles.plotStyle = 'bars';
end
        
function scatterBtnCB(hObj,event)
  % check whether plot style is already boxWhisker
  if strcmp(handles.plotStyle,'1Dscat')
    return
  end
  % reverse border on each
  changeBorder(handles.boxWhiskJava,0);
  changeBorder(handles.barJava,0);
  changeBorder(handles.scatterJava,1);
  % change plot style
  handles.plotStyle = '1Dscat';
end

function closeCB(hObj,event)
  close(gcf);  
end

function runCB(hObj,event)
%   if ~checkControls()
%     return
%   end
  
  kitLog('Starting intensity analysis');
  
  % check what data type has been selected
  dataTypes = {'intrameas','jobset'};
  type = handles.dataType{2}.Value+1;
  type = dataTypes{type};
  
  % get channel information and labels
  contchan = handles.controlChNum;
  analchans = zeros(1,3);
  for iChan = 1:3
    analchans(iChan) = handles.analyseCh{iChan}.Value;
    chanlabel{iChan}  = handles.chlabels{iChan}.String;
  end
  analchans = find(analchans);
  chanlabel = chanlabel(analchans);
  
  % get all jobs ready for analysis
  nConds = length(handles.tabP.Children);
  jobs = repmat({[]},1,nConds);
  
  % predesignate analysis structure if using intra-measurements
  if strcmp(type,'intrameas')
    norms = struct('cellwise',[],'spotwise',[]);
    analysis = struct('channels',[1 2],...
                      'raw',[],...
                      'norm',norms,...
                      'bg',[],...
                      'stats',[]);
    analysis = repmat({analysis},1,nConds);
  end
  
  % loop over conditions
  for iCond = 1:nConds    
    % find which, if any, jobsets are selected
    selData = handles.jobsets{iCond}.Value;
    if isempty(selData)
      errorbox('Must first select data in every condition before executing analysis.');
      return
    end
    datdir = handles.jobsetdir{iCond}.String;
    
    % get condition information
    condlabel{iCond} = handles.condlabel{iCond}.String;
    normcond(iCond) = handles.normConds{iCond}.Value;
    
    % get intensities per condition, using wither data type
    switch type
      
      case 'jobset'
        
        % loop over selected jobsets for this condition
        for iData = selData
          jsname = handles.jobsets{iCond}.String{iData};
          jobset = kitLoadJobset(fullfile(datdir,jsname));
          tempjobs = kitLoadAllJobs(jobset);
          jobs{iCond} = [jobs{iCond} tempjobs];
        end
        kitLog('Getting intensities: %s',condlabel{iCond});
    
        % get intensity measurements
        analysis{iCond} = intensityMeasurements(jobs{iCond},...
            'controlChan',contchan,'analysisChans',analchans);
          
      case 'intrameas'
        
        for iData = selData
          
          % get intrameasurement data from file
          imname = handles.jobsets{iCond}.String{iData};
          imdata = load(fullfile(datdir,imname));
          fnames = fieldnames(imdata);
          imdata = getfield(imdata,fnames{1});
              
          % collate the data per experiment
          filtIdx = ~isnan(imdata.microscope.depthFilter.delta.threeD.all(:));
          tempRaw = [imdata.intensity.mean.inner(:) imdata.intensity.mean.outer(:)];
          tempRaw = tempRaw(filtIdx,:); %filter out data not included in delta measurements
          analysis{iCond}.raw = [analysis{iCond}.raw; tempRaw];
          tempSpot = tempRaw;
          tempSpot(:,1) = tempSpot(:,1)./tempRaw(:,contchan);
          tempSpot(:,2) = tempSpot(:,2)./tempRaw(:,contchan);
          analysis{iCond}.norm.spotwise = [analysis{iCond}.norm.spotwise; tempSpot];
          analysis{iCond}.norm.cellwise = [analysis{iCond}.norm.cellwise; nan(size(tempSpot))];
          tempBg = [imdata.intensity.bg.inner(:) imdata.intensity.bg.outer(:)];
          analysis{iCond}.bg = [analysis{iCond}.bg; tempBg];
            
        end
              
      otherwise
        errorbox('Ensure only jobsets or intra-measurements structures are selected.');  
        return
    end
      
  end
  kitLog('Intensities compiled across all conditions');
  
  % get condition information for plotting
  normcond = find(normcond);
  plotstyle = handles.plotStyle;
  closefigs = ~handles.closeFigs.Value;
  % get saving filenames etc
  filename = handles.filename.String;
  savepath = handles.savedir;
  
  % create filename
  GUIfile = ['GUIoptions_' filename];
  % print the figure to file
  print(fullfile(savepath,GUIfile),'-depsc');
  
  % produce the figures
  intensityFigures(analysis,'channelLabels',chanlabel,'conditions',condlabel,...
      'normalise',normcond,'save',1,'savePath',savepath,'filename',filename,...
      'closeFigs',closefigs,'plotType',plotstyle);
  
  kitLog('Spot detection analysis complete')
  
end

%% Other functions.

function populateJSbox(condidx)
  jobsetDir = handles.jobsetdir{condidx}.String;
  if ~isempty(jobsetDir)
    % Find jobset files.
    waitmsg = 'Finding jobsets...';
    hwait = waitbar(0,waitmsg);
    if handles.allfiles{condidx}.Value
      jobsets = kitFindFiles(jobsetDir,'.mat');
    else
      jobsets = kitFindFiles(jobsetDir,'.mat');
      rmv = contains(jobsets,[filesep 'kittracking']);
      jobsets(rmv) = [];
    end
    % Strip search directory from filenames.
    for i=1:length(jobsets)
      jobsets{i} = strrep(jobsets{i},[jobsetDir filesep],'');
    end
    set(handles.jobsets{condidx},'String',jobsets,'Value',1:length(jobsets));
    waitbar(1,hwait,waitmsg);
    close(hwait);
  end
end

function changeBorder(javaH,add)
  % set up new border details
  rndCrnrs = true;
  col = java.awt.Color(1,0,0);
  if add
    thickness = 1;
  else
    thickness = 0;
  end
  newBorder = javax.swing.border.LineBorder(col,thickness,rndCrnrs);
  % give new border details to the java handle
  set(javaH,'Border',newBorder);
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

end % kitIntensityAnalysis
