function jobset=kitGUI(jobset)
% KITGUI Display GUI to setup and run tracking.
%
% Copyright (c) 2015 Jonathan W. Armond

% Download BioFormats, if required.
kitDownloadBioFormats();

if nargin<1 || isempty(jobset)
  jobset = kitDefaultOptions();
end

% Upgrade jobset, if required.
if ~isfield(jobset,'jobsetVersion') || ...
    jobset.jobsetVersion < kitVersion(2)
  jobset = kitJobset(jobset);
end

if ~isfield(jobset,'ROI')
  jobset.ROI = [];
end

coordSystemValues = {'Plate','Image','Centre of mass'}; % in GUI
coordSystemValuesJS = {'plate','image','com'}; % in jobset
spotDetectValues = {'Histogram','Adaptive','Wavelet','Manual','None'};
spotDetectValuesJS = {'histcut','adaptive','wavelet','manual','none'};
spotRefineValues = {'Centroid','MMF','None'};
spotRefineValuesJS = {'centroid','gaussian','none'};

% Setup GUI.
handles = createControls();
updateControls(jobset);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls()
  colwidth = [55 42 35];
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 sum(colwidth)+19 44]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiT ' kitVersion(1)];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  figpos = hs.fig.Position;
  figw = figpos(3);
  figh = figpos(4);
  headfont = 16;
  medfont = 14;
  smallfont = 12;

  w=25; h=8;
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[figpos(3)-w-3 1 w h]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kitlogo.png'),pos([4 3])));

  %% ROI selection
  x = 2.5;
  w = colwidth(1);
  toplabely = figh-2;
  % Movies
  hs.openBtn = button(hs.fig,'Open existing...',[x toplabely 20 2],@openExistingCB);
  hs.movies = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[x 0.52*figh w 0.34*figh],'Max',inf,'Min',0,'FontSize',smallfont);
  hs.movies.String = jobset.movieFiles;
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.movies);
  jScrollPane.setHorizontalScrollBarPolicy(32);

  hs.movieDirectory = editbox(hs.fig,'',[x 0.91*figh w 1.7]);
  hs.movieDirectory.Enable = 'inactive';
  hs.movieDirectory.String = jobset.movieDirectory;
  hs.selectDirectory = button(hs.fig,'Select directory',[w-15 0.865*figh 17.5 2],@selectDirectoryCB);
  hs.labelAvail = label(hs.fig,'Available movies:',[x 0.871*figh 20 1.5]);

  % ROIs
  hs.ROIs = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[x 3 w 18.15],'Callback',@roisCB);
  hs.ROIs.Min = 1;
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.ROIs);
  jScrollPane.setHorizontalScrollBarPolicy(32);

  label(hs.fig,'ROIs:',[x 21.2 20 1.5]);
  hs.addROI = button(hs.fig,'Add ROI',[x 1 13 2],@addROICB);
  hs.viewROI = button(hs.fig,'View ROI',[x+14.5 1 13 2],@viewROICB);
  hs.deleteROI = button(hs.fig,'Delete ROI',[x+14.5*2 1 13 2],@deleteROICB);

  %% Tracking setup
  x = colwidth(1) + 3;
  w = colwidth(2);
  t = label(hs.fig,'Tracking setup',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  % Coordinate system
  label(hs.fig,'Coordinate system',[x 40 20 1.5]);
  hs.coordSys = popup(hs.fig,coordSystemValues,[x+23 40.1 20 1.5]);
  label(hs.fig,'Coordinate system channel',[x 38 30 1.5]);
  hs.coordSysCh = editbox(hs.fig,[],[x+35 38 6 1.5]);

  % Channel modes
  b = 32;
  h = 5.5;
  for i=1:3
    p = uipanel(hs.fig,'Units','characters','Position',[x b-(i-1)*h w h],'FontSize',12,'Title',['Channel ' num2str(i)]);
    hs.spotMode{i} = popup(p,spotDetectValues,[22 2.5 0.45*w 1.5],@spotModeCB);
    label(p,'Spot detection',[1 2.5 20 1.5]);

    hs.refineMode{i} = popup(p,spotRefineValues,[22 0.5 0.45*w 1.5],@refineModeCB);
    label(p,'Spot refinement',[1 0.5 20 1.5]);
  end

  b = b-2*h;
  % Tasks
  % tasks = {'Spot finding','Plane fitting','Tracking','Sister grouping','Intensity measurement'};
  % h = 1.75;
  % panelh = h*(1+length(tasks));
  % p = uipanel(hs.fig,'Units','characters','Position',[x b-panelh w panelh],'FontSize',12,'Title','Tasks');
  % for i=1:length(tasks)
  %    hs.tasks{i} = checkbox(p,tasks{i},[1 panelh-h*(i+1) 30 h]);
  %    hs.tasks{i}.Value = hs.tasks{i}.Max;
  % end
  % b = b-panelh;

  %% Execution
  y = b-h;
  labelw = 0.5*w;
  t = label(hs.fig,'Execution',[x y labelw 1.5],14);
  t.FontWeight = 'bold';
  h = 2;
  y = y-h;
  label(hs.fig,'Jobset name',[x y labelw h],12);
  hs.filename = editbox(hs.fig,'jobset',[x+w-(w-labelw) y (w-labelw) h]);
  btnw = 0.5*w;
  bx = x + w - btnw;
  y = y-h;
  hs.validateMetadata = button(hs.fig,'Validate metadata',[bx y btnw h],@validateCB);
  y = y-h;
  hs.save = button(hs.fig,'Save',[bx y btnw h],@saveCB);
  y = y-h;
  hs.execute = button(hs.fig,'Execute',[bx y btnw h],@executeCB);
  y = y-h;
  hs.parallel = checkbox(hs.fig,'Execute in parallel',[x y w h]);
  if ~license('test','Distrib_Computing_Toolbox') || verLessThan('distcomp','6.5')
    hs.parallel.Enable = 'off';
  end

  %% Options
  x = sum(colwidth(1:2)) + 4;
  w = 35;
  t = label(hs.fig,'Tracking options',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  h = 1.5;
  lh = 1.5*h; % large height
  y = toplabely-h;
  hs.autoRadii = checkbox(hs.fig,'Calculate search radii from dt',[x y w h],@autoRadiiCB,10);
  y = y-h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  t = label(hs.fig,'Frame dt',[x y labelw h],10);
  hs.autoRadiidt = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Est. avg. disp. of spots (um/s)',[x y labelw h],10);
  % Assume mean absolute displacment of sisters is about 0.06 μm/s as default.
  hs.autoRadiiAvgDisp = editbox(hs.fig,num2str(0.06),[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Min search radius (um)',[x y labelw h],10);
  hs.minSearchRadius = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Max search radius (um)',[x y labelw h],10);
  hs.maxSearchRadius = editbox(hs.fig,[],[editx y editw h],10);

  y = y-h;
  hs.useSisterAlignment = checkbox(hs.fig,'Use sister alignment',[x y w h],@useSisterAlignmentCB,10);
  y = y-h;
  % Adjust text box pos for multiple lines.
  t = label(hs.fig,'Max angle between sisters and plate normal (deg)',[x y-h/2 labelw lh],10);
  hs.maxSisterAlignmentAngle = editbox(hs.fig,[],[editx y editw h],10);
  y = y-lh;
  t = label(hs.fig,'Max average distance between sisters (um)',[x y-h/2 labelw lh],10);
  hs.maxSisterDist = editbox(hs.fig,[],[editx y editw h],10);
  y = y-lh;
  t = label(hs.fig,'Min overlap between sister tracks',[x y-h/2 labelw lh],10);
  hs.minSisterTrackOverlap = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Number of manual detect frames',[x y labelw h],10);
  hs.nManualDetectFrames = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Min spots per frame',[x y labelw h],10);
  hs.minSpotsPerFrame = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Weight for spot count',[x y labelw h],10);
  hs.adaptiveLambda = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.mmfAddSpots = checkbox(hs.fig,'Resolve sub-resolution spots',[x y w h],'',10);
  y = y-h;
  t = label(hs.fig,'Max MMF time per frame (min)',[x y labelw h],10);
  hs.maxMmfTime = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.alphaAText = label(hs.fig,'Weight for intensity restriction:',[x y labelw h],10);
  hs.alphaAchText{1} = label(hs.fig,'Ch.1',[x+labelw-3 y-0.15 4 h],10);
  hs.alphaA{1} = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.alphaAchText{2} = label(hs.fig,'Ch.2',[x+labelw-3 y-0.15 4 h],10);
  hs.alphaA{2} = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.alphaAchText{3} = label(hs.fig,'Ch.3',[x+labelw-3 y-0.15 4 h],10);
  hs.alphaA{3} = editbox(hs.fig,[],[editx y editw h],10);
%   y = y-h;
%   hs.deconvolve = checkbox(hs.fig,'Deconvolve movies',[x y w h],@deconvolveCB,10);
%   y = y-h;
%   hs.psfBtn = button(hs.fig,'Select PSF',[x y labelw/2 h],@psfBtnCB,10);
%   hs.psfFile = editbox(hs.fig,'',[x+labelw/2+2 y colwidth(3)-labelw/2-2 h],10);

  movegui(hs.fig,'center');
end

% Update control status based on contents of jobset.
function updateControls(jobset)
  hs = handles;
  opts = jobset.options;

  if isfield(jobset,'movieDirectory')
    handles.movieDirectory.String = jobset.movieDirectory;
  end
  if isfield(jobset,'filename')
    [~,file] = fileparts(jobset.filename);
    hs.filename.String = file;
  end
  hs.coordSys.Value = mapStrings(opts.coordSystem,coordSystemValuesJS);
  hs.coordSysCh.String = num2str(opts.coordSystemChannel);
  for i=1:3
    hs.spotMode{i}.Value = mapStrings(opts.spotMode{i},spotDetectValuesJS);
    hs.refineMode{i}.Value = mapStrings(opts.coordMode{i},spotRefineValuesJS);
  end
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
  hs.nManualDetectFrames.String = num2str(opts.manualDetect.numFrames);
  hs.minSpotsPerFrame.String = num2str(opts.minSpotsPerFrame);
  hs.adaptiveLambda.String = num2str(opts.adaptiveLambda);
  hs.mmfAddSpots.Value = opts.mmfAddSpots;
  hs.maxMmfTime.String = num2str(opts.maxMmfTime);
  for iChan=1:3;
    hs.alphaA{iChan}.String = num2str(opts.alphaA(iChan));
  end
  if isfield(jobset,'psfFile')
    hs.psfFile.String = jobset.psfFile;
  end
%   if opts.deconvolve
%     hs.deconvolve.Value = hs.deconvolve.Max; % On
%     hs.psfBtn.Enable = 'on';
%     hs.psfFile.Enable = 'on';
%   else
%     hs.deconvolve.Value = hs.deconvolve.Min; % Off.
%     hs.psfBtn.Enable = 'off';
%     hs.psfFile.Enable = 'off';
%   end

  populateMovieBox();
  populateROIBox();
  spotModeCB();
  refineModeCB();
  autoRadiiCB();
end

% Check controls for consistent input.
function tf=checkControls()
  hs = handles;
  v = str2double(hs.coordSysCh.String);
  if ~isfinite(v) || v<1 || v>3
    errorbox('Invalid channel number for coordinate system channel. Should be a number between 1 and 3');
    tf = false;
    return
  end

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

function roisCB(hObj,event)
  if isempty(handles.ROIs.Value)
    handles.ROIs.Value = 1;
  end
end

function openExistingCB(hObj,event)
  if ~isempty(get(handles.ROIs,'String'))
    r = questdlg('Selecting existing jobset will clear existing ROIs. Select?','Warning','Yes','No','No');
    if strcmp(r,'No')
      return
    end
  end

  [filename,pathname] = uigetfile('*.mat','Select existing jobset');
  if ~isempty(filename)
    filename = fullfile(pathname,filename);
    try
      jobset = kitLoadJobset(filename);
      if isempty(jobset) || ~isfield(jobset,'kit') || ~jobset.kit
        error('Jobset file corrupt');
      end
      % Upgrade jobset, if required.
      if ~isfield(jobset,'jobsetVersion') || ...
          jobset.jobsetVersion < kitVersion(2)
        jobset = kitJobset(jobset);
      end

      if ~isfield(jobset,'ROI')
        jobset.ROI = [];
      end

      updateControls(jobset);
      h=msgbox('Successfully loaded jobset','Success','Help','modal');
      uiwait(h);
    catch me
      errorbox(sprintf('Error loading jobset %s: %s',filename,me.message));
    end
  end
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
    ROI = [];
  end
end

function addROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  if isempty(v)
    errorbox('Must select movies first to add ROIs');
    return
  end
  for i=1:length(v)
    movieFileName = fullfile(movieDir,movieFiles{v(i)});
    kitLog('Opening movie file %i of %i: %s',i,length(v),movieFileName);
    [crop,cropSize] = kitCropMovie(movieFileName);
    if ~isempty(crop)
      for j=1:size(crop,1)
        r = length(jobset.ROI) + 1;
        jobset.ROI(r).movie = handles.movies.String{v(i)};
        jobset.ROI(r).crop = crop(j,:);
        jobset.ROI(r).cropSize = cropSize(j,:);
      end
    end
  end
  populateROIBox();
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

function viewROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    movieDir = handles.movieDirectory.String;
    kitMovieProj(fullfile(movieDir,jobset.ROI(v).movie),[],jobset.ROI(v).crop);
  end
end

function spotModeCB(hObj,event)
  if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotDetectValues),'Adaptive'),handles.spotMode))
    handles.adaptiveLambda.Enable = 'on';
  else
    handles.adaptiveLambda.Enable = 'off';
  end
  if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotDetectValues),'Manual'),handles.spotMode))
    handles.nManualDetectFrames.Enable = 'on';
  else
    handles.nManualDetectFrames.Enable = 'off';
  end
end

function refineModeCB(hObj,event)
  if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotRefineValues),'MMF'),handles.refineMode))
    handles.mmfAddSpots.Enable = 'on';
    handles.maxMmfTime.Enable = 'on';
    for iChan=1:3
      if strcmp(mapStrings(handles.refineMode{iChan}.Value,spotRefineValues),'MMF');
        handles.alphaA{iChan}.Enable = 'on';
      else
        handles.alphaA{iChan}.Enable = 'off';
      end
    end
  else
    handles.mmfAddSpots.Enable = 'off';
    handles.maxMmfTime.Enable = 'off';
    for iChan=1:3;
      handles.alphaA{iChan}.Enable = 'off';
    end
  end
end

function autoRadiiCB(hObj,event)
  if handles.autoRadii.Value
    handles.autoRadiidt.Enable = 'on';
    handles.autoRadiiAvgDisp.Enable = 'on';
    handles.minSearchRadius.Enable = 'off';
    handles.maxSearchRadius.Enable = 'off';
    r = computeSearchRadii(str2double(handles.autoRadiidt.String),str2double(handles.autoRadiiAvgDisp.String));
    handles.minSearchRadius.String = num2str(r(1));
    handles.maxSearchRadius.String = num2str(r(2));
  else
    handles.autoRadiidt.Enable = 'off';
    handles.autoRadiiAvgDisp.Enable = 'off';
    handles.minSearchRadius.Enable = 'on';
    handles.maxSearchRadius.Enable = 'on';
  end
end

function useSisterAlignmentCB(hObj,event)
  if handles.useSisterAlignment.Value
    handles.maxSisterAlignmentAngle.Enable = 'on';
  else
    handles.maxSisterAlignmentAngle.Enable = 'off';
  end
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

function executeCB(hObj,event)
  if ~checkControls()
    return
  end
  updateJobset();
  kitSaveJobset(jobset);
  % Ask which tasks to run.
  taskStrs = {'Spot finding','Plane fitting','Tracking','Sister grouping','Intensity measurement'};
  [tasks,ok] = listdlg('ListString',taskStrs,'InitialValue',1:length(taskStrs),'PromptString','Select tasks to execute','Name','Select tasks...');
  if ok
    % Map to task numbers and add defaults.
    taskMap = [1 2 3 4 8];
    tasks = [taskMap(tasks) 5:7];

    if handles.parallel.Value
      execmode = 'batch';
    else
      execmode = 'serial';
    end
    progh = waitbar(0,sprintf('Tracking progress (%d/%d)',0,length(jobset.ROI)));
    kitRunJobs(jobset,'callback',@trackProgress,'tasks',tasks,'exec',execmode);
    delete(progh);
    uiresume(gcf);
  end

  function trackProgress(idx)
    waitbar(idx/length(jobset.ROI),progh,sprintf('Tracking progress (%d/%d)',idx,length(jobset.ROI)));
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

% function deconvolveCB(hObj,event)
%   if handles.deconvolve.Value
%     handles.psfBtn.Enable = 'on';
%     handles.psfFile.Enable = 'on';
%   else
%     handles.psfBtn.Enable = 'off';
%     handles.psfFile.Enable = 'off';
%   end
% end

% function psfBtnCB(hObj,event)
%   [file,path] = uigetfile('*.mat','Select PSF file');
%   if isequal(file,0)
%     return
%   end
%   handles.psfFile.String = file;
%   file = fullfile(path,file);
%   data = load(file);
%   % Assume PSF is only variable.
%   f = fieldnames(data);
%   if length(f) > 1
%     h=msgbox(sprintf('Multiple variables in MAT-file. Using .%s for PSF.',f{1}),'Warning','Warning','modal');
%     uiwait(h);
%   end
%   jobset.psf = data.(f{1});
%   jobset.psfFile = file;
% end

function populateMovieBox()
  movieDir = handles.movieDirectory.String;
  if isempty(movieDir)
    handles.movies.String = [];
  else
    % Find movie files.
    movieFiles = kitFindFiles(movieDir, kitSupportedFormats(),0,1);
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
  if (handles.ROIs.Value <= 0 && length(handles.ROIs.String)>0) || handles.ROIs.Value > length(jobset.ROI)
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
  jobset.filename = fullfile(jobset.movieDirectory,[handles.filename.String '.mat']);

  opts = jobset.options;
  opts.coordSystem = mapStrings(handles.coordSys.Value,coordSystemValuesJS);
  opts.coordSystemChannel = str2double(handles.coordSysCh.String);
  for i=1:3
    opts.spotMode{i} = spotDetectValuesJS{handles.spotMode{i}.Value};
    opts.coordMode{i} = spotRefineValuesJS{handles.refineMode{i}.Value};
  end
  if handles.autoRadii.Value
    opts.autoRadiidt = str2double(handles.autoRadiidt.String);
    opts.autoRadiiAvgDisp = str2double(handles.autoRadiiAvgDisp.String);
    r = computeSearchRadii(opts.autoRadiidt,opts.autoRadiiAvgDisp);
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
  opts.manualDetect.numFrames = str2double(handles.nManualDetectFrames.String);
  opts.minSpotsPerFrame = str2double(handles.minSpotsPerFrame.String);
  opts.adaptiveLambda = str2double(handles.adaptiveLambda.String);
  opts.mmfAddSpots = handles.mmfAddSpots.Value;
  opts.maxMmfTime = str2double(handles.maxMmfTime.String);
  for iChan=1:3;
    opts.alphaA(iChan) = str2double(handles.alphaA{iChan}.String);
  end
%   opts.deconvolve = handles.deconvolve.Value;
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

end % kitGUI
