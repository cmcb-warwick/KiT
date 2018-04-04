function jobset=kitSetupJob(jobset)
% KITSETUPJOB Display a GUI to set up a job.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

% Download BioFormats, if required.
kitDownloadBioFormats();

if nargin<1 || isempty(jobset)
  jobset = kitDefaultOptions();
  handles.mode = chooseMode;
else
  handles.mode = jobset.options.jobProcess;
end
jobset.options.jobProcess = handles.mode;

% Upgrade jobset, if required.
if ~isfield(jobset,'jobsetVersion') || ...
    jobset.jobsetVersion < kitVersion(2)
  jobset = kitJobset(jobset);
end

if ~isfield(jobset,'ROI')
  jobset.ROI = [];
end

coordSystemValues = {'Plane fit','Image moments','Centre of mass','None'};
coordSystemValuesJS = {'plate','image','com','none'};
spotDetectValues = {'Histogram','Adaptive','Wavelet','Manual'};
spotDetectValuesJS = {'histcut','adaptive','wavelet','manual'};
spotRefineValues = {'Centroid','MMF','None'};
spotRefineValuesJS = {'centroid','gaussian','none'};
maskValues = {'Circle','Semi-circle','Cone'};
maskValuesJS = {'circle','semicirc','cone'};

% Setup GUI.
handles = createControls();
updateControls(jobset);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls()
  
  % Give mode to the new hs structure.
  hs.mode = handles.mode;
    
  % Set up figure.  
  colw = [55 50];
  figw = sum(colw);
  figh = 49;
  figpos = [100 35 figw figh];
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',figpos);
  set(0,'DefaultFigureWindowStyle','normal');
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  switch hs.mode
      case 'zandt'
        hs.fig.Name = 'Set up new tracking job';
      case 'zonly'
        hs.fig.Name = 'Set up new single timepoint job';
      case 'chrshift'
        hs.fig.Name = 'Set up new chromatic shift job';
  end
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  % define some font sizes
  largefont = 16;
  medfont = 14;
  smallfont = 12;
  tinyfont = 10;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 0.5; %small horizontal shift
  toplabely = figh; %top-most point

  % Movie selection sub-title.
  labx = dx; y = toplabely-lh;
  labw = colw(1)-2*dx;
  t = label(hs.fig,'1. Movie and ROI selection',[labx y labw h],largefont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Directory selection.
  btnx = dx;
  btnw = 15; btnh = 2;
  hs.selectDirectory = button(hs.fig,'Select directory',[btnx y btnw btnh],@selectDirectoryCB);
  labx = dx+(btnw+ddx); labw = colw(1)-(labx+dx);
  hs.movieDirectory = editbox(hs.fig,'',[labx y labw h]);
  hs.movieDirectory.Enable = 'inactive';
  hs.movieDirectory.String = jobset.movieDirectory;
  
  % Movies panel
  labx = dx; y = y-lh;
  labw = colw(1)-2*dx;
  hs.labelAvail = label(hs.fig,'All available movies:',[labx y labw h]);
  panw = colw(1)-2*dx; panh = 12*h;
  panx = dx; y = y-panh;
  hs.movies = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[panx y panw panh],'Max',inf,'Min',0,'FontSize',smallfont);
  hs.movies.String = jobset.movieFiles;
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.movies);
  jScrollPane.setHorizontalScrollBarPolicy(32);
  
  % ROI label and buttons.
  y = y-lh;
  labx = dx; labw = colw(1)-2*dx;
  label(hs.fig,'Selected movies:',[labx y labw h]);
  btnx = dx; y = y-lh;
  btnw = (colw(1)-(btnx+dx+3*ddx))/4; btnh = 2;
  hs.fullROI = button(hs.fig,'Add',[btnx y btnw btnh],@addROICB,smallfont);
  btnx = btnx+(btnw+ddx);
  hs.cropROI = button(hs.fig,'Add + crop',[btnx y btnw btnh],@cropROICB,smallfont);
  btnx = btnx+(btnw+ddx);
  hs.deleteROI = button(hs.fig,'Delete',[btnx y btnw btnh],@deleteROICB,smallfont);
  btnx = btnx+(btnw+ddx);
  hs.viewROI = button(hs.fig,'View',[btnx y btnw btnh],@viewROICB,smallfont);
  
  % ROI panel.
  panw = colw(1)-2*dx; panh = 12*h;
  panx = dx; y = y-(panh+ddx);
  hs.ROIs = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[panx y panw panh],'Callback',@roisCB);
  hs.ROIs.Min = 1;
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.ROIs);
  jScrollPane.setHorizontalScrollBarPolicy(32);

  % Channel setup sub-title.
  labx = colw(1)+dx; y = toplabely-lh;
  labw = colw(2)-2*dx;
  t = label(hs.fig,'2. Process setup',[labx y labw h],largefont);
  t.FontWeight = 'bold';
  
  % Set some consistent label and edit thicknesses.
  labx = colw(1)+dx; labw = 25;
  editx = labx+(labw+ddx); editw = colw(2)-(labw+ddx+2*dx);
  
  if ~strcmp(hs.mode,'chrshift')
    % Coordinate system
    y = y-lh;
    hs.coordSysText = label(hs.fig,'Coordinate system',[labx y labw h]);
    hs.coordSys = popup(hs.fig,coordSystemValues,[editx y editw h],@neighbourChCB);
  end
  
  % Detection channel information.
  y = y-lh;
  label(hs.fig,'Detect spots in channel...',[labx y labw h]);
  radx = editx; radw = editw/4;
  for i=1:3
    hs.spotDetectCh{i} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String',num2str(i),...
        'Position',[radx y radw h],'Callback',@spotDetectChCB,'FontSize',tinyfont);
    radx = radx+radw;
  end
  y = y-h;
  hs.spotDetectCh{1}.Value = 1;
  hs.spotDetectChNum = 1;
  
  % Detection and refinement methods.
  label(hs.fig,'Spot detection',[labx y labw h]);
  hs.detectMode = popup(hs.fig,spotDetectValues,[editx y editw h],@detectModeCB);
  y = y-h;
  label(hs.fig,'Spot refinement',[labx y labw h]);
  hs.refineMode = popup(hs.fig,spotRefineValues,[editx y editw h],@refineModeCB);
  
  % Neighbour channel information.
  y = y-lh;
  hs.neighbourChText = label(hs.fig,'Use channel 1 to detect spots in channels...',[labx y colw(2)-2*dx h]);
  y = y-h;
  radx = editx; radw = editw/4;
  for i=1:3
    hs.neighbourCh{i} = checkbox(hs.fig,num2str(i),[radx y radw h],@neighbourChCB,tinyfont);
    radx = radx+radw;
  end
  hs.neighbourCh{2}.Value = 1;
  hs.neighbourCh{1}.Enable = 'off';

  % Options sub-title.
  labx = colw(1)+dx; y = y-lh;
  labw = colw(2)-2*dx;
  t = label(hs.fig,'3. Options',[labx y labw h],largefont);
  t.FontWeight = 'bold';
  
  % Create options tabs.
  panw = colw(2)-2*dx; panh = 13.5*h;
  panx = colw(1)+dx; pany = y-panh;
  hs.tabPanel = uipanel(hs.fig,'Units','characters',...
      'Position',[panx pany panw panh],'FontSize',smallfont,'BorderType','none');
  hs.tabP = uitabgroup('Parent', hs.tabPanel);
  
  % Set some standard positions and tallies for tabs.
  tabID = 0;
  tabtoplabely = panh-lh;
  w = panw-2*dx;
  labx = dx; labw = 30;
  editx = labx+(labw+ddx); editw = w-editx;
  
  % Detection tab.
  tabID = tabID+1;
  hs.tabs{tabID} = uitab('Parent', hs.tabP, 'Title', 'Detection');
  taby = tabtoplabely;
  t = label(hs.tabs{tabID},'Primary spot detection',[labx taby labw h],smallfont);
  t.FontWeight = 'Bold';
  taby = taby-h;
  hs.minSpotsText = label(hs.tabs{tabID},'Min spots per frame',[labx taby labw h],tinyfont);
  hs.minSpots = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  taby = taby-h;
  hs.maxSpotsText = label(hs.tabs{tabID},'Max spots per frame',[labx taby labw h],tinyfont);
  hs.maxSpots = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  taby = taby-h;
  hs.manualFrameSpaceText = label(hs.tabs{tabID},'Manual detection frame spacing',[labx taby labw h],tinyfont);
  hs.manualFrameSpace = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  taby = taby-lh;
  t = label(hs.tabs{tabID},'Secondary spot detection',[labx taby labw h],smallfont);
  t.FontWeight = 'Bold';
  taby = taby-h;
  hs.neighbourMaskShapeText = label(hs.tabs{tabID},'Mask shape',[labx taby labw h],tinyfont);
  hs.neighbourMaskShape = popup(hs.tabs{tabID},maskValues,[editx-editw taby editw*2 h],@neighbourChCB,tinyfont);
  taby = taby-h;
  hs.neighbourMaskRadiusText = label(hs.tabs{tabID},'Mask radius (um)',[labx taby labw h],tinyfont);
  hs.neighbourMaskRadius = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  if ~strcmp(hs.mode,'chrshift')
      taby = taby-h;
      hs.chOrientText = label(hs.tabs{tabID},'Channel orientation',[labx taby w h],tinyfont);
      taby = taby-h;
      labw = (w-2*dx)/8.4; labx = 2*dx;
      logow = labw*16/25; logoh = labw*4/15;
      hs.chOrientInner = label(hs.tabs{tabID},'inner KT',[labx taby-h/4 labw lh],tinyfont);
      labx = labx+labw;
      for i=1:3
        hs.chOrient{i} = label(hs.tabs{tabID},num2str(i),[labx taby labw h],smallfont);
        hs.chOrient{i}.HorizontalAlignment = 'center';
        labx = labx+labw;
        if i==3
          continue
        end
        hs.chSwap{i} = button(hs.tabs{tabID},num2str(i),[labx taby+h/4 logow logoh],@chSwapCB);
        hs.chSwap{i}.ForegroundColor = [0.94 0.94 0.94];
        pos = getpixelposition(hs.chSwap{i});
        [I,~]=imread('private/channelswap.png','BackgroundColor',[0.94 0.94 0.94]);
        set(hs.chSwap{i},'cdata',imresize(I,pos([4 3])),'HorizontalAlignment','center');
        labx = labx+logow;
      end
      hs.chOrientOuter = label(hs.tabs{tabID},'outer KT',[labx taby-h/4 labw lh],tinyfont);
      labx = dx; labw = 30; % reset some values
  end
  
  % Mixture model fitting tab.
  tabID = tabID+1;
  hs.tabs{tabID} = uitab('Parent', hs.tabP, 'Title', 'MMF');
  taby = tabtoplabely;
  hs.mmfAddSpots = checkbox(hs.tabs{tabID},'Resolve sub-resolution spots',[labx taby w h],'',tinyfont);
  taby = taby-h;
  hs.maxMmfTimeText = label(hs.tabs{tabID},'Max MMF time per frame (s)',[labx taby labw h],tinyfont);
  hs.maxMmfTime = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  taby = taby-lh;
  hs.alphaAText = label(hs.tabs{tabID},'Weight for intensity restriction in channel...',[labx taby labw h],tinyfont);
  for i=1:3
    hs.alphaAtext{i} = label(hs.tabs{tabID},num2str(i),[editx-(2+ddx) taby 2 h],tinyfont);
    hs.alphaAtext{i}.HorizontalAlignment = 'right';
    hs.alphaA{i} = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
    taby = taby-h;
  end
  
  if strcmp(hs.mode,'zandt')
      % Tracking tab.
      tabID = tabID+1;
      hs.tabs{tabID} = uitab('Parent', hs.tabP, 'Title', 'Tracking');
      taby = tabtoplabely;
      hs.autoRadii = checkbox(hs.tabs{tabID},'Calculate search radii from dt',[labx taby w h],@autoRadiiCB,tinyfont);
      taby = taby-h;
      hs.autoRadiidtText = label(hs.tabs{tabID},'Time lapse (s)',[labx taby labw h],tinyfont);
      hs.autoRadiidt = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
      taby = taby-h;
      hs.autoRadiiAvgDispText = label(hs.tabs{tabID},'Est. avg. kinetochore displacement (um/min)',[labx taby labw h],tinyfont);
      % Assume mean absolute displacment of sisters is about 3.6 um/min as default.
      hs.autoRadiiAvgDisp = editbox(hs.tabs{tabID},num2str(3.6),[editx taby editw h],tinyfont);
      taby = taby-h;
      hs.minSearchRadiusText = label(hs.tabs{tabID},'Min search radius (um)',[labx taby labw h],tinyfont);
      hs.minSearchRadius = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
      taby = taby-h;
      hs.maxSearchRadiusText = label(hs.tabs{tabID},'Max search radius (um)',[labx taby labw h],tinyfont);
      hs.maxSearchRadius = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
      taby = taby-h;
      hs.useSisterAlignment = checkbox(hs.tabs{tabID},'Use sister alignment',[labx taby w h],@useSisterAlignmentCB,tinyfont);
      taby = taby-h;
      hs.maxSisterAlignmentAngleText = label(hs.tabs{tabID},'Max angle between sisters and plate normal (deg)',[labx taby-h/2 labw lh],tinyfont);
      hs.maxSisterAlignmentAngle = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
      taby = taby-lh;
      hs.maxSisterDistText = label(hs.tabs{tabID},'Max average distance between sisters (um)',[labx taby-h/2 labw lh],tinyfont);
      hs.maxSisterDist = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
      taby = taby-h;
      hs.minSisterTrackOverlapText = label(hs.tabs{tabID},'Min overlap between sister tracks',[labx taby labw h],tinyfont);
      hs.minSisterTrackOverlap = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  end
  
  if ~strcmp(hs.mode,'chrshift')
      % Chromatic shift tab.
      tabID = tabID+1;
      hs.tabs{tabID} = uitab('Parent', hs.tabP, 'Title', 'Chr. shift');
      taby = tabtoplabely;
      hs.chromaticShift = checkbox(hs.tabs{tabID},'Provide chromatic shift correction',[labx taby w h],@chromaticShiftCB,tinyfont);
      cspanh = 6*h;
      taby = taby-cspanh;
      hs.csPanel = uipanel(hs.tabs{tabID},'Units','characters','Position',[labx taby w-dx cspanh],'FontSize',tinyfont,'Title','Chromatic shift jobsets');
      p = hs.csPanel;
      taby = taby-lh;
      % Create jobset panel within tab.
      labw = 6; btnw = 15.5; editw = 2;
      labx = dx;
      btnx = labx+(labw+dx);
      editx1 = btnx+(btnw+dx);
      editx2 = editx1+(2*editw);
      cspany = cspanh-(lh+h);
      hs.csVectText = label(p,'Channel vector',[labx cspany labw lh],tinyfont);
      hs.csVectText.FontWeight = 'bold';
      hs.csJobsetText = label(p,'Jobset',[btnx cspany labw h],tinyfont);
      hs.csJobsetText.FontWeight = 'bold';
      hs.csOrderText = label(p,'Channel order',[editx1 cspany labw lh],tinyfont);
      hs.csOrderText.FontWeight = 'bold';
      cspany = cspany-h;
      for iChan = 1:(3-1)
        hs.csVect{iChan} = label(p,sprintf('1  ->  %i',iChan+1),[labx cspany labw h],tinyfont);
        hs.csArrow{iChan} = label(p,'->',[(editx1+editx2)/2 cspany editw h],tinyfont);
        hs.csArrow{iChan}.HorizontalAlignment = 'center';
        hs.csJobset{iChan} = button(p,'-',[btnx cspany btnw h],@csOrderCB,tinyfont);
        hs.csOrder{iChan,1} = editbox(p,[],[editx1 cspany editw h],tinyfont);
        hs.csOrder{iChan,2} = editbox(p,[],[editx2 cspany editw h],tinyfont);
        cspany = cspany-h;
      end
      labw = 30; editw = w-editx;
      hs.csFilter = checkbox(hs.tabs{tabID},'Filter chromatic shift spots',[labx taby w h],@csFilterCB,tinyfont);
      taby = taby-h;
      hs.csAmplitudeText = label(hs.tabs{tabID},'Min spot intensity (% of max)',[labx taby labw h],tinyfont);
      hs.csAmplitude = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
      taby = taby-h;
      hs.csDistanceText = label(hs.tabs{tabID},'Min spot separation (um)',[labx taby labw h],tinyfont);
      hs.csDistance = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  else
      hs.chromaticShift.Value = 0;
  end
  
  if ~strcmp(hs.mode,'chrshift')
      % Intensity tab.
      tabID = tabID+1;
      hs.tabs{tabID} = uitab('Parent', hs.tabP, 'Title', 'Intensity');
      taby = tabtoplabely;
      label(hs.tabs{tabID},'Measure in channels...',[labx taby labw h],tinyfont);
      labw = 17.5;
      radx = labx+(labw+ddx); radw = (w-radx-3*ddx)/4;
      for i=1:3
        hs.intensityExecute{i} = checkbox(hs.tabs{tabID},num2str(i),[radx taby radw h],@intensityOptionsCB,tinyfont);
        radx = radx+(radw+ddx);
      end
      labw = 30;
      taby = taby-h;
      hs.intensityMaskShapeText = label(hs.tabs{tabID},'Mask shape',[labx taby labw h],tinyfont);
      hs.intensityMaskShape = popup(hs.tabs{tabID},maskValues,[editx-editw taby editw*2 h],[],tinyfont);
      taby = taby-h;
      hs.intensityMaskRadiusText = label(hs.tabs{tabID},'Mask radius (um)',[labx taby labw h],tinyfont);
      hs.intensityMaskRadius = editbox(hs.tabs{tabID},[],[editx taby editw h],tinyfont);
  end
  
  % Tasks tab.
  tabID = tabID+1;
  hs.tabs{tabID} = uitab('Parent', hs.tabP, 'Title', 'Tasks');
  tasks = {'Primary spot detection','Coordinate system fitting','Spot tracking','Sister pairing','Secondary spot detection','Intensity measurement'};
  switch hs.mode
      case 'zonly'
          tasks = tasks([1 2 5 6]);
      case 'chrshift'
          tasks = tasks([1 5]);
  end
  hs.taskList = tasks;
  taby = tabtoplabely;
  for i=1:length(tasks)
     hs.tasks{i} = checkbox(hs.tabs{tabID},tasks{i},[dx taby labw h]);
     hs.tasks{i}.Value = hs.tasks{i}.Max;
     taby = taby-h;
  end
  
  % Reassign w.
  w = colw(2)-2*dx;  
  y = y-(panh+h);

  % Final checks title.
  labx = colw(1)+dx; labw = w;
  t = label(hs.fig,'4. Final checks',[labx y labw h],largefont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Jobset name.
  labw = w; 
  t = label(hs.fig,'Filename',[labx y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-h; labw = 15;
  hs.filenameTxt = label(hs.fig,['kitjobset_' datestr(now,'yymmdd') '_'],[labx+ddx y labw h],smallfont);
  hs.filename = editbox(hs.fig,'filename',[labx+(labw+ddx) y w-(labw+ddx) h]);
%   y = y-h;
%   tickw = w;
%   hs.groupOutput = checkbox(hs.fig,'Group output files into separate folder',[labx y tickw h],[],tinyfont);
  y = y-lh;
  btnw = 17.5; btnx = labx+(w/2-btnw/2); btnh = 2;
  hs.validateMetadata = button(hs.fig,'Validate metadata',[btnx y btnw btnh],@validateCB);
  y = y-btnh;
  hs.save = button(hs.fig,'Save jobset',[btnx y btnw btnh],@saveCB);
  
  % Cancel button.
  btnw = 12; btnh = 2;
  x = figw-(btnw+dx); y = 1;
  hs.cancel = button(hs.fig,'Cancel',[x y btnw btnh],@cancelCB);
  
  movegui(hs.fig,'center');
  
end

% Update control status based on contents of jobset.
function updateControls(jobset)
  hs = handles;
  opts = jobset.options;
  
  % check kitdetection grouping option
%   hs.groupOutput.Value = opts.groupOutput;

  if isfield(jobset,'movieDirectory')
    handles.movieDirectory.String = jobset.movieDirectory;
  end
  if isfield(jobset,'filename')
    [~,file] = fileparts(jobset.filename);
    file = file(18:end);
    hs.filename.String = file;
  end
  if ~strcmp(hs.mode,'chrshift')
    hs.coordSys.Value = mapStrings(opts.coordSystem,coordSystemValuesJS);
  end
  hs.spotDetectChNum = opts.coordSystemChannel;
  for i=1:3
    hs.spotDetectCh{i}.Value = (i==hs.spotDetectChNum);
    hs.neighbourCh{i}.Value = strcmp(opts.spotMode{i},'neighbour');
  end
  hs.detectMode.Value = mapStrings(opts.spotMode{opts.coordSystemChannel},spotDetectValuesJS);
  hs.refineMode.Value = mapStrings(opts.coordMode{opts.coordSystemChannel},spotRefineValuesJS);
  if strcmp(hs.mode,'zandt')
    hs.autoRadiidt.String = num2str(opts.autoRadiidt);
    if isempty(opts.autoRadiidt)
      hs.autoRadii.Value = hs.autoRadii.Min; % Off
      hs.autoRadiidt.Enable = 'off';
      hs.autoRadiiAvgDisp.Enable = 'off';
      hs.minSearchRadius.Enable = 'on';
      hs.maxSearchRadius.Enable = 'on';
    else
      hs.autoRadii.Value = hs.autoRadii.Max; % On
      hs.autoRadiidt.Enable = 'on';
      hs.autoRadiiAvgDisp.Enable = 'on';
      hs.minSearchRadius.Enable = 'off';
      hs.maxSearchRadius.Enable = 'off';
    end
    hs.minSearchRadius.String = num2str(opts.minSearchRadius(1));
    hs.maxSearchRadius.String = num2str(opts.maxSearchRadius(1));
    hs.useSisterAlignment.Value = opts.useSisterAlignment;
    hs.maxSisterAlignmentAngle.String = num2str(opts.maxSisterAlignmentAngle);
    if ~hs.useSisterAlignment.Value
      hs.maxSisterAlignmentAngle.Enable = 'off';
    end
    hs.maxSisterDist.String = num2str(opts.maxSisterSeparation);
    hs.minSisterTrackOverlap.String = num2str(opts.minSisterTrackOverlap);
    autoRadiiCB();
  end
  
  hs.minSpots.String = num2str(opts.minSpotsPerFrame);
  hs.maxSpots.String = num2str(opts.maxSpotsPerFrame);
  hs.manualFrameSpace.String = num2str(opts.manualDetect.frameSpacing);
  hs.mmfAddSpots.Value = opts.mmf.addSpots;
  hs.maxMmfTime.String = num2str(opts.mmf.maxMmfTime);
  for iChan=1:3
    hs.alphaA{iChan}.String = num2str(opts.mmf.alphaA(iChan));
  end
  
  if ~strcmp(hs.mode,'chrshift')
    
    hs.chromaticShift.Value = any(~cellfun('isempty',opts.chrShift.jobset(:)));
    hs.csMinSpots.String = num2str(opts.chrShift.minSpots);
    % list channels and select in order
    chans = setdiff(1:3,opts.coordSystemChannel,'stable');
    for chanID = 1:(3-1)
      iChan = chans(chanID);
      if ~isempty(opts.chrShift.jobset{1,iChan})
        hs.csJobset{chanID}.String = opts.chrShift.jobset{1,iChan};
        hs.csOrder{chanID,1}.Enable = 'on';
        hs.csOrder{chanID,1}.String = num2str(opts.chrShift.chanOrder{1,iChan}(1));
        hs.csOrder{chanID,2}.Enable = 'on';
        hs.csOrder{chanID,2}.String = num2str(opts.chrShift.chanOrder{1,iChan}(2));
      else
        hs.csJobset{chanID}.String = '-';
        hs.csOrder{chanID,1}.String = '';
        hs.csOrder{chanID,1}.Enable = 'off';
        hs.csOrder{chanID,2}.String = '';
        hs.csOrder{chanID,2}.Enable = 'off';
      end
    end
    hs.csFilter.Value = opts.chrShift.filtering;
    hs.csAmplitude.String = num2str(opts.chrShift.intensityFilter*100);
    hs.csDistance.String = num2str(opts.chrShift.neighbourFilter);
    chromaticShiftCB();
    csFilterCB();
    
    for iChan=1:3
      hs.intensityExecute{iChan}.Value = opts.intensity.execute(iChan);
    end
    hs.intensityMaskShape.Value = mapStrings(opts.intensity.maskShape,maskValuesJS);
    hs.intensityMaskRadius.String = num2str(opts.intensity.maskRadius);
    intensityOptionsCB();
    
    for iChan=1:3
      hs.chOrient{iChan}.String = num2str(opts.neighbourSpots.channelOrientation(iChan));
    end
    hs.neighbourMaskShape.Value = mapStrings(opts.neighbourSpots.maskShape,maskValuesJS);
    
  else
    hs.neighbourMaskShape.Value = mapStrings('circle',maskValuesJS);
  end
  hs.neighbourMaskRadius.String = num2str(opts.neighbourSpots.maskRadius);
  
  if isfield(jobset,'metadata')
      hs.validateMetadata.String = 'Re-validate...';
  end
  
  handles = hs;
  
  populateMovieBox();
  populateROIBox();
  spotDetectChCB();
  detectModeCB();
  refineModeCB();
  neighbourChCB();
  chromaticShiftCB();
  
end

% Check controls for consistent input.
function tf=checkControls()
  
  hs = handles;

  if strcmp(hs.mode,'zandt')
    v = str2double(hs.autoRadiidt.String);
    if (hs.autoRadii.Value == hs.autoRadii.Max) && (~isfinite(v) || v<=0)
      errorbox('Invalid value for frame dt. Should be a positive number.')
      tf = false;
      return
    end

    v = str2double(hs.autoRadiiAvgDisp.String);
    if (hs.autoRadii.Value == hs.autoRadii.Max) && (~isfinite(v) || v<0)
      errorbox('Invalid value for avg. disp. of spots. Should be a positive number.')
      tf = false;
      return
    end

    v1 = str2double(hs.minSearchRadius.String);
    v2 = str2double(hs.maxSearchRadius.String);
    if (hs.autoRadii.Value ~= hs.autoRadii.Max) && (~isfinite(v1) || v1 > v2 || v1 < 0)
      errorbox('Invalid value min search radius. Should be a positive number less than max search radius.')
      tf = false;
      return
    end
    if (hs.autoRadii.Value ~= hs.autoRadii.Max) && (~isfinite(v2) || v2 < v1 || v2 < 0)
      errorbox('Invalid value max search radius. Should be a positive number less than min search radius.')
      tf = false;
      return
    end

    v = str2double(hs.maxSisterAlignmentAngle.String);
    if (hs.useSisterAlignment.Value == hs.useSisterAlignment.Max) && (~isfinite(v) || v < 0 || v > 180)
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
    
  end

  v = str2double(hs.minSpots.String);
  if ~isfinite(v) || v < 0
    errorbox('Invalid value for min spots per frame. Should be a positive number.')
    tf = false;
    return
  end
  v(2) = str2double(hs.maxSpots.String);
  if ~isfinite(v(2)) || v(2) < 0
    errorbox('Invalid value for maximum spots per frame. Should be a positive number.')
    tf = false;
    return
  elseif diff(v)<=0
    errorbox('Invalid values for spots per frame. Maximum number should be larger than the minimum.')
    tf = false;
    return
  end

  v = str2double(hs.maxMmfTime.String);
  if (hs.mmfAddSpots.Value == hs.mmfAddSpots.Max) && (~isfinite(v) || v < 0)
    errorbox('Invalid value for min spots per frame. Should be a positive number or zero.')
    tf = false;
    return
  end

  if isempty(hs.filename.String)
    errorbox('Jobset name is a required field.');
    tf = false;
    return
  end

  tf = true;
end

function selectDirectoryCB(hObj,event)
  if ~isempty(get(handles.ROIs,'String'))
    r = questdlg('Selecting movie directory will clear existing ROIs. Select?','Warning','Yes','No','No');
    if strcmp(r,'No')
      return
    end
  end

  dirName = uigetdir([], 'Select directory tree containing movies');
  if ~isempty(dirName)
    set(handles.movieDirectory, 'String', dirName);
    populateMovieBox();
    set(handles.ROIs,'String',[]);
  end
end

function roisCB(hObj,event)
  if isempty(handles.ROIs.Value)
    handles.ROIs.Value = 1;
  end
end

function cropROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  if isempty(v)
    errorbox('Must select movies first to add ROIs');
    return
  end
  % Start waitbar.
  waitmsg = 'Finding jobsets...';
  hwait = waitbar(0,waitmsg);
  % Loop over selected movies.
  for i=1:length(v)
    [crop,cropSize] = kitCropMovie(fullfile(movieDir,movieFiles{v(i)}));
    if ~isempty(crop)
      for j=1:size(crop,1)
        r = length(jobset.ROI) + 1;
        jobset.ROI(r).movie = handles.movies.String{v(i)};
        jobset.ROI(r).crop = crop(j,:);
        jobset.ROI(r).cropSize = cropSize(j,:);
      end
    end
    waitbar(i/length(v),waitmsg);
  end
  populateROIBox();
  close(hwait);
end

function deleteROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    r = questdlg('Delete selected ROI?','Warning','Yes','No','No');
    if strcmp(r,'Yes')
      jobset.ROI(v) = [];
    end
  end
  populateROIBox();
end

function addROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  if isempty(v)
    errorbox('Must select movies first to skip ROIs');
    return
  end
  % Start waitbar.
  waitmsg = 'Cropping movies...';
  hwait = waitbar(0,waitmsg);
  % Loop over selected movies.
  for i=1:length(v)
    [md,~]=kitOpenMovie(fullfile(movieDir,movieFiles{v(i)}));
    crop = [1 1 md.frameSize(1:2)];
    cropSize = md.frameSize(1:3);
    r = length(jobset.ROI) + 1;
    jobset.ROI(r).movie = handles.movies.String{v(i)};
    jobset.ROI(r).crop = crop;
    jobset.ROI(r).cropSize = cropSize;
    waitbar(i/length(v),hwait,waitmsg);
  end
  populateROIBox();
  close(hwait);
end

function viewROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    movieDir = handles.movieDirectory.String;
    kitMovieProj(fullfile(movieDir,jobset.ROI(v).movie),[],jobset.ROI(v).crop);
  end
end

function spotDetectChCB(hObj,event)
  if exist('hObj','var')
    chan = str2double(hObj.String);
    handles.spotDetectChNum = chan;
  else
    chan = handles.spotDetectChNum;
  end
  handles.spotDetectCh{chan}.Value = 1;
  handles.neighbourCh{chan}.Value = 0;
  handles.neighbourCh{chan}.Enable = 'off';
  for notChan = setdiff(1:3,chan)
    handles.spotDetectCh{notChan}.Value = 0;
    handles.neighbourCh{notChan}.Enable = 'on';
  end
  handles.neighbourChText.String(13) = num2str(chan);
  neighbourChCB();
end

function detectModeCB(hObj,event)
  if strcmp(mapStrings(handles.detectMode.Value,spotDetectValues),'Histogram')
    handles.minSpots.Enable = 'on';
    handles.minSpotsText.Enable = 'on';
    handles.maxSpots.Enable = 'on';
    handles.maxSpotsText.Enable = 'on';
  else
    handles.minSpots.Enable = 'off';
    handles.minSpotsText.Enable = 'off';
    handles.maxSpots.Enable = 'off';
    handles.maxSpotsText.Enable = 'off';
  end
  if strcmp(mapStrings(handles.detectMode.Value,spotDetectValues),'Manual')
    handles.manualFrameSpaceText.Enable = 'on';
    handles.manualFrameSpace.Enable = 'on';
  else
    handles.manualFrameSpaceText.Enable = 'off';
    handles.manualFrameSpace.Enable = 'off';
  end
  neighbourChCB();
end

function refineModeCB(hObj,event)
  if strcmp(mapStrings(handles.refineMode.Value,spotRefineValues),'MMF')
    handles.mmfAddSpots.Enable = 'on';
    handles.maxMmfTimeText.Enable = 'on';
    handles.maxMmfTime.Enable = 'on';
    handles.alphaAText.Enable = 'on';
    for iChan=1:3
      if handles.neighbourCh{iChan}.Value || iChan == handles.spotDetectChNum
        handles.alphaAtext{iChan}.Enable = 'on';
        handles.alphaA{iChan}.Enable = 'on';
      else
        handles.alphaAtext{iChan}.Enable = 'off';
        handles.alphaA{iChan}.Enable = 'off';
      end
    end
  else
    handles.mmfAddSpots.Enable = 'off';
    handles.maxMmfTimeText.Enable = 'off';
    handles.maxMmfTime.Enable = 'off';
    handles.alphaAText.Enable = 'off';
    for iChan=1:3
      handles.alphaAtext{iChan}.Enable = 'off';
      handles.alphaA{iChan}.Enable = 'off';
    end
  end
end

function neighbourChCB(hObj,event)
  neighChans = cellfun(@(x) x.Value,handles.neighbourCh);
  if sum(neighChans)>0
    handles.neighbourMaskRadiusText.Enable = 'on';
    handles.neighbourMaskRadius.Enable = 'on';
    if ~strcmp(handles.mode,'chrshift')
        handles.neighbourMaskShapeText.Enable = 'on';
        handles.neighbourMaskShape.Enable = 'on';
        handles.chOrientText.Enable = 'on';
        for i=1:3
          handles.chOrient{i}.Enable = 'on';
          if i~=3
            handles.chSwap{i}.Visible = 'on';
          end
        end
        handles.chOrientInner.Enable = 'on';
        handles.chOrientOuter.Enable = 'on';
    else
        handles.neighbourMaskShapeText.Enable = 'off';
        handles.neighbourMaskShape.Enable = 'off';
    end
  else
    handles.neighbourMaskShapeText.Enable = 'off';
    handles.neighbourMaskShape.Enable = 'off';
    handles.neighbourMaskRadiusText.Enable = 'off';
    handles.neighbourMaskRadius.Enable = 'off';
    if ~strcmp(handles.mode,'chrshift')
        handles.chOrientText.Enable = 'off';
        for i=1:3
          handles.chOrient{i}.Enable = 'off';
          if i~=3
            handles.chSwap{i}.Visible = 'off';
          end
        end
        handles.chOrientInner.Enable = 'off';
        handles.chOrientOuter.Enable = 'off';
    end
  end
  
  if ~strcmp(handles.mode,'chrshift') && any(strcmp(mapStrings(handles.coordSys.Value,coordSystemValues),{'Centre of mass','None'}))
    handles.neighbourMaskShapeText.Enable = 'off';
    handles.neighbourMaskShape.Value = 1;
    handles.neighbourMaskShape.Enable = 'off';
    handles.chOrientText.Enable = 'off';
    for i=1:3
      handles.chOrient{i}.Enable = 'off';
      if i~=3
        handles.chSwap{i}.Visible = 'off';
      end
    end
    handles.chOrientInner.Enable = 'off';
    handles.chOrientOuter.Enable = 'off';
  end
  refineModeCB();
  taskOptionsCB();
  
end

function autoRadiiCB(hObj,event)
  if handles.autoRadii.Value
    handles.autoRadiidtText.Enable = 'on';
    handles.autoRadiidt.Enable = 'on';
    handles.autoRadiiAvgDispText.Enable = 'on';
    handles.autoRadiiAvgDisp.Enable = 'on';
    handles.minSearchRadiusText.Enable = 'off';
    handles.minSearchRadius.Enable = 'off';
    handles.maxSearchRadiusText.Enable = 'off';
    handles.maxSearchRadius.Enable = 'off';
    r = computeSearchRadii(str2double(handles.autoRadiidt.String),str2double(handles.autoRadiiAvgDisp.String)/60);
    handles.minSearchRadius.String = num2str(r(1));
    handles.maxSearchRadius.String = num2str(r(2));
  else
    handles.autoRadiidtText.Enable = 'off';
    handles.autoRadiidt.Enable = 'off';
    handles.autoRadiiAvgDispText.Enable = 'off';
    handles.autoRadiiAvgDisp.Enable = 'off';
    handles.minSearchRadiusText.Enable = 'on';
    handles.minSearchRadius.Enable = 'on';
    handles.maxSearchRadiusText.Enable = 'on';
    handles.maxSearchRadius.Enable = 'on';
  end
end

function useSisterAlignmentCB(hObj,event)
  if handles.useSisterAlignment.Value
    handles.maxSisterAlignmentAngleText.Enable = 'on';
    handles.maxSisterAlignmentAngle.Enable = 'on';
  else
    handles.maxSisterAlignmentAngleText.Enable = 'off';  
    handles.maxSisterAlignmentAngle.Enable = 'off';
  end
end

function chSwapCB(hObj,event)
  refPos = str2double(hObj.String);
  pos{1} = handles.chOrient{refPos}.String;
  pos{2} = handles.chOrient{refPos+1}.String;
  handles.chOrient{refPos}.String = pos{2};
  handles.chOrient{refPos+1}.String = pos{1};
end

function chromaticShiftCB(hObj,event)
  if handles.chromaticShift.Value
    handles.csMinSpotsText.Enable = 'on';
    handles.csMinSpots.Enable = 'on';
    handles.csPanel.ForegroundColor = [0 0 0];
    handles.csVectText.Enable = 'on';
    handles.csJobsetText.Enable = 'on';
    handles.csOrderText.Enable = 'on';
    for iChan = 1:(3-1)
      handles.csVect{iChan}.Enable = 'on';
      handles.csArrow{iChan}.Enable = 'on';
      handles.csJobset{iChan}.Enable = 'on';
      if any(strcmp(handles.csJobset{iChan}.String,{'-','Unknown source'}))
        handles.csOrder{iChan,1}.Enable = 'on';
        handles.csOrder{iChan,2}.Enable = 'on';
        handles.csOrder{iChan,1}.String = '1';
        handles.csOrder{iChan,2}.String = num2str(iChan+1);
      else
        handles.csOrder{iChan,1}.Enable = 'off';
        handles.csOrder{iChan,2}.Enable = 'off';
        handles.csOrder{iChan,1}.String = '';
        handles.csOrder{iChan,2}.String = '';
      end
    end
    handles.csFilter.Enable = 'on';
  else
    handles.csMinSpotsText.Enable = 'off';
    handles.csMinSpots.Enable = 'off';
    handles.csPanel.ForegroundColor = [0.5 0.5 0.5];
    handles.csVectText.Enable = 'off';
    handles.csJobsetText.Enable = 'off';
    handles.csOrderText.Enable = 'off';
    for iChan = 1:(3-1)
      handles.csVect{iChan}.Enable = 'off';
      handles.csArrow{iChan}.Enable = 'off';
      handles.csJobset{iChan}.Enable = 'off';
      handles.csOrder{iChan,1}.Enable = 'off';
      handles.csOrder{iChan,2}.Enable = 'off';
      handles.csOrder{iChan,1}.String = '';
      handles.csOrder{iChan,2}.String = '';
    end
    handles.csFilter.Value = 0;
    csFilterCB();
    handles.csFilter.Enable = 'off';
  end
end

function csFilterCB(hObj,event)
  if handles.csFilter.Value
    handles.csAmplitudeText.Enable = 'on';
    handles.csAmplitude.Enable = 'on';
    handles.csDistanceText.Enable = 'on';
    handles.csDistance.Enable = 'on';
  else
    handles.csAmplitudeText.Enable = 'off';
    handles.csAmplitude.Enable = 'off';
    handles.csDistanceText.Enable = 'off';
    handles.csDistance.Enable = 'off';
  end
end

function csOrderCB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  for i=1:(3-1)
    testPos(i) = handles.csJobset{i}.Position(2);
  end
  chan = find(testPos == hObj.Position(2));
  handles.csJobset{chan}.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{1,chan+1} = file;
  handles.csOrder{chan,1}.Enable = 'on';
  handles.csOrder{chan,1}.String = jobset.options.chrShift.chanOrder{1,chan+1}(1);
  handles.csOrder{chan,2}.Enable = 'on';
  handles.csOrder{chan,2}.String = jobset.options.chrShift.chanOrder{1,chan+1}(2);
end

function intensityOptionsCB(hObj,event)
  
  if any(cellfun(@(x) x.Value,handles.intensityExecute))
      handles.intensityMaskShapeText.Enable = 'on';
      handles.intensityMaskShape.Enable = 'on';
      handles.intensityMaskRadiusText.Enable = 'on';
      handles.intensityMaskRadius.Enable = 'on';
  else
      handles.intensityMaskShapeText.Enable = 'off';
      handles.intensityMaskShape.Enable = 'off';
      handles.intensityMaskRadiusText.Enable = 'off';
      handles.intensityMaskRadius.Enable = 'off';
  end
  taskOptionsCB();
    
end

function taskOptionsCB(hObj,event)
  
  % Neighbour spot detection.
  taskID = find(cellfun(@(x) strcmp(x,'Secondary spot detection'),handles.taskList));
  if any(cellfun(@(x) x.Value,handles.neighbourCh))
      handles.tasks{taskID}.Value = 1;
      handles.tasks{taskID}.Enable = 'on';
  else
      handles.tasks{taskID}.Value = 0;
      handles.tasks{taskID}.Enable = 'off';
  end
  
  if ~strcmp(handles.mode,'chrshift')
    % Intensity measurement.
    taskID = find(cellfun(@(x) strcmp(x,'Intensity measurement'),handles.taskList));
    if any(cellfun(@(x) x.Value,handles.intensityExecute))
      handles.tasks{taskID}.Value = 1;
      handles.tasks{taskID}.Enable = 'on';
    else
      handles.tasks{taskID}.Value = 0;
      handles.tasks{taskID}.Enable = 'off';
    end
  end
    
end

function saveCB(hObj,event)
  if ~checkControls()
    return
  end
  updateJobset();
  kitSaveJobset(jobset);
  uiresume(gcf);
end

function validateCB(hObj,event)
  updateJobset();
  % check that ROIs are already provided
  if ~isfield(jobset,'ROI')
    errorbox('Must provide ROIs prior to validating their metadata.')
  else
    % loop over each movie
    for iMov = 1:length(jobset.ROI)
      [jobset,applyAll] = kitValidateMetadata(jobset,iMov);
      % skip showing remaining movies, just save metadata for all
      if applyAll
        for jMov = iMov+1:length(jobset.ROI)
          jobset.metadata{jMov} = jobset.metadata{iMov};
        end
        handles.validateMetadata.String = 'Re-validate...';
        break
      end
    end
    
  end
  % show that movies have been validated
  handles.validateMetadata.String = 'Re-validate...';
  
end

function cancelCB(hObj,event)
  jobset.cancel = 1;
  uiresume(gcf);
end

%% Other functions

function populateMovieBox()
  movieDir = handles.movieDirectory.String;
  if isempty(movieDir)
    handles.movies.String = [];
  else
    % Find movie files.
    movieFiles = kitFindFiles(movieDir, kitSupportedFormats(),1,0,1);
    % Strip search directory from filenames.
    for i=1:length(movieFiles)
      movieFiles{i} = strrep(movieFiles{i},[movieDir filesep],'');
    end
    set(handles.movies,'String',movieFiles,'Value',1:length(movieFiles));
  end
end

function populateROIBox()
  handles.ROIs.String=[];
  maxMovLen = 32;
  movieFiles = handles.movies.String;
  handles.ROIs.Value = 0; % Keep MATLAB quiet about invalid selection.
  if ~isempty(movieFiles)
    for i=1:length(jobset.ROI)
      handles.ROIs.String{i} = [strshorten(jobset.ROI(i).movie,maxMovLen) ' [' ...
                          num2str(round(jobset.ROI(i).crop),'%d ') ']'];
    end
  end
  if (handles.ROIs.Value <= 0 && ~isempty(handles.ROIs.String)) || handles.ROIs.Value > length(jobset.ROI)
    handles.ROIs.Value = 1;
  end
end

function ret=mapStrings(inp,vals)
  if ischar(inp)
    % Map string to index.
    ret = find(cellfun(@(x) ~isempty(x),strfind(vals,inp)));
    assert(~isempty(ret));
  else
    % Map index to index.
    assert(inp > 0 && inp <= length(vals));
    ret = vals{inp};
  end
end

function updateJobset()
  jobset.movieDirectory = handles.movieDirectory.String;
  jobset.movieFiles = handles.movies.String;
  
  filename = [handles.filenameTxt.String handles.filename.String '.mat'];
  jobset.filename = fullfile(jobset.movieDirectory,filename);

  opts = jobset.options;
  
  % Process options.
  if ~strcmp(handles.mode,'chrshift')
    opts.coordSystem = mapStrings(handles.coordSys.Value,coordSystemValuesJS);
  else
    opts.coordSystem = 'com';
  end
  opts.coordSystemChannel = handles.spotDetectChNum;
  opts.spotMode{handles.spotDetectChNum} = spotDetectValuesJS{handles.detectMode.Value};
  opts.coordMode{handles.spotDetectChNum} = spotRefineValuesJS{handles.refineMode.Value};
  for i=setdiff(1:3,handles.spotDetectChNum)
    if handles.neighbourCh{i}.Value
      opts.spotMode{i} = 'neighbour';
      opts.coordMode{i} = spotRefineValuesJS{handles.refineMode.Value};
    else
      opts.spotMode{i} = 'none';
      opts.coordMode{i} = 'none';
    end
  end
  if strcmp(handles.mode,'zandt')
    % Tracking options.
    if handles.autoRadii.Value
      opts.autoRadiidt = str2double(handles.autoRadiidt.String);
%       opts.autoRadiiAvgDisp = str2double(handles.autoRadiiAvgDisp.String)/60;
      r = computeSearchRadii(opts.autoRadiidt,str2double(handles.autoRadiiAvgDisp)/60);
    else
      opts.autoRadiidt = [];
      r = zeros(2,1);
      r(1) = str2double(handles.minSearchRadius.String);
      r(2) = str2double(handles.maxSearchRadius.String);
    end
    r = computeUnalignedLaggingRadii(r);
    opts.minSearchRadius = r(1,:);
    opts.maxSearchRadius = r(2,:); % in um
    opts.useSisterAlignment = handles.useSisterAlignment.Value;
    opts.maxSisterAlignmentAngle = str2double(handles.maxSisterAlignmentAngle.String);
    opts.maxSisterSeparation = str2double(handles.maxSisterDist.String);
    opts.minSisterTrackOverlap = str2double(handles.minSisterTrackOverlap.String);
  end
  % Primary spot detection options.
  opts.minSpotsPerFrame = str2double(handles.minSpots.String);
  opts.maxSpotsPerFrame = str2double(handles.maxSpots.String);
  opts.manualDetect.frameSpacing = str2double(handles.manualFrameSpace.String);
  % MMF options.
  mmf = opts.mmf;
  mmf.addSpots = handles.mmfAddSpots.Value;
  mmf.maxMmfTime = str2double(handles.maxMmfTime.String);
  for iChan=1:3
    mmf.alphaA(iChan) = str2double(handles.alphaA{iChan}.String);
  end
  opts.mmf = mmf;
  % Chromatic shift options.
  if handles.chromaticShift.Value
    chrShift = opts.chrShift;
    chrShift.minSpots = str2double(handles.csMinSpots.String);
    chans = setdiff(1:3,handles.spotDetectChNum);
    for iChan = 1:(3-1)
      jChan = chans(iChan);
      chrShift.chanOrder{opts.coordSystemChannel,jChan}(1) = str2double(handles.csOrder{iChan,1}.String);
      chrShift.chanOrder{opts.coordSystemChannel,jChan}(2) = str2double(handles.csOrder{iChan,2}.String);
    end
    if handles.csFilter.Value
      chrShift.filtering = 1;
      chrShift.intensityFilter = str2double(handles.csAmplitude.String)/100;
      chrShift.neighbourFilter = str2double(handles.csDistance.String);
    end
    result = getChromaticShiftResults(chrShift);
    chrShift.result = result;
    opts.chrShift = chrShift;
  end
  % Neighbour spot detection options.
  neighbourSpots = opts.neighbourSpots;
  neighbourSpots.maskShape = mapStrings(handles.neighbourMaskShape.Value,maskValuesJS);
  neighbourSpots.maskRadius = str2double(handles.neighbourMaskRadius.String);
  if ~strcmp(handles.mode,'chrshift')
    for iChan=1:3
      neighbourSpots.channelOrientation(iChan) = str2double(handles.chOrient{iChan}.String);
    end
  end
  opts.neighbourSpots = neighbourSpots;
  if ~strcmp(handles.mode,'chrshift')
    % Intensity options.
    intensity = opts.intensity;
    for iChan=1:3
      intensity.execute(iChan) = handles.intensityExecute{iChan}.Value;
    end
    intensity.maskShape = mapStrings(handles.intensityMaskShape.Value,maskValuesJS);
    intensity.maskRadius = str2double(handles.intensityMaskRadius.String);
    intensity.photobleachCorrect = intensity.photobleachCorrect*strcmp(handles.mode,'zandt');
    opts.intensity = intensity;
  end
  
%   opts.groupOutput = handles.groupOutput.Value;
  
  jobset.options = opts;
end

function r=computeSearchRadii(dt,avgDisp)
  r = zeros(2,1);
  r(2) = 6*avgDisp*dt;
  r(1) = 0.1*avgDisp*dt;
end

function r=computeUnalignedLaggingRadii(r)
% Assume unaligned move 3x faster, lagging same speed.
  r = [r 3*r r];
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

function cellResult = getChromaticShiftResults(chrShift)
  cellResult = chrShift.result;
  for i = 1:3
    for j = 1:3
      jS = chrShift.jobset{i,j};
	  if i>=j || isempty(jS) || strcmp(jS,'Unknown source'); continue; end
      jS = kitLoadJobset(jS);
      neighFilt = chrShift.neighbourFilter;
      intFilt = chrShift.intensityFilter/100;
      mS = kitLoadAllJobs(jS);
      if chrShift.filtering
        mS = chrsFilterSpots(mS,'revert',1, ...
          'neighbourFilter',neighFilt,'intensityFilter',intFilt, ...
          'referenceChan',handles.spotDetectChNum);
      end
      [result,~] = chrsCalculateChromaticShift(mS,[i j],...
          'filtered',1);
      cellResult{i,j} = result; cellResult{j,i} = result.*[-1 -1 -1 1 1 1];
    end
  end       
end

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

end % kitGUI
