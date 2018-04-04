function jobset=kitChangeOptions(jobset,mode)
% KITCHANGEOPTIONS Displays GUI to allow user to change JOBSET options.
%
% For debugging, use jobset=kitChangeOptions(jobset,'debug').
%
% Copyright (c) 2017 C. A. Smith

if nargin<1 || isempty(jobset)
  openExistingCB();
end

neighbourMaskValues = {'Circle','Semi-circle','Cone'};
neighbourMaskValuesJS = {'circle','semicirc','cone'};
directionMethodValues = {'Voting','Step-wise'};
directionMethodValuesJS = {'voting','stepwise'};
debugGroupSisterValues = {'-','for 4 frames','for all tracks'};
debugMMFclustersValues = {'-','without pausing','with a pause','and stop'};
debugMMFcandsValues = {'-','without pausing','with a pause'};
debugMMFfinalValues = {'-','without pausing','with a pause'};

% get options framework
opts = jobset.options;
tabsToShow = ones(1,5);
    
% check whether or not debug has been requested
tabsToShow(1) = ~(nargin<2 || isempty(mode) || ~strcmp(mode,'debug'));

% check whether or not chromatic shift is appropriate here
tabsToShow(2) = (sum(cellfun(@(x) ~strcmp(x,'none'),opts.spotMode))>1 ...
      && ~strcmp(opts.jobProcess,'chrshift'));
  
% check whether MMF is required for any channel
tabsToShow(4) = any(cellfun(@(x) strcmp(x,'gaussian'),opts.coordMode));

% check whether tracking is required here
tabsToShow(5) = strcmp(opts.jobProcess,'zandt');

% Setup GUI.
handles = createControls();
updateControls(opts);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls()
  
  colwidth = 55;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 colwidth+5 40]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiT ' kitVersion(1) ' - Update options'];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  figpos = hs.fig.Position;
  figw = figpos(3);
  figh = figpos(4);
  headfont = 16;
  medfont = 14;
  smallfont = 12;
  tinyfont = 10;
  
  % Create tabs
  hs.tabs = uitabgroup('Parent', hs.fig);
  if tabsToShow(1)
    hs.debugTab = uitab('Parent', hs.tabs, 'Title', 'Debugging');
  end
  if tabsToShow(2)
    hs.chrShiftTab = uitab('Parent', hs.tabs, 'Title', 'Chromatic shift');
  end
  if tabsToShow(3)
    hs.spotDetectTab = uitab('Parent', hs.tabs, 'Title', 'Spot detection');
  end
  if tabsToShow(4)
    hs.mmfTab = uitab('Parent', hs.tabs, 'Title', 'MMF');
  end
  if tabsToShow(5)
    hs.trackingTab = uitab('Parent', hs.tabs, 'Title', 'Tracking');
  end
  
  % Set some standard positions and distances
  h = 1.5; %height
  lh = 1.5*h; %large height
  x = 2.5;
  toplabely = figh-2*x; %top-most point
  
  % Editable jobset name at top of the screen
  hs.jobsetName = editbox(hs.fig,[],[x toplabely+0.75 figw-2*x h],medfont);
  
  % Save and 'Save as...' buttons
  btnw = 10;
  btnx = figw-2.5-btnw;
  hs.update = button(hs.fig,'Save',[btnx 3 btnw h],@saveCB);
  hs.newJobset = button(hs.fig,'Save as...',[btnx 1.25 btnw h],@newJobsetCB);
  
  % KiT logo
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[btnx 4.75 btnw 3.2]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kitlogo.png'),pos([4 3])));
  
  
  %% Chromatic shift options, tab 1
  
  if tabsToShow(2)
  
  % tab-specific positions and distances
  colwidth = 35;
  x = x+10;
  w = colwidth-5;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  y = toplabely-lh;
  
  % whether or not to provide chromatic shift
  hs.chromaticShift = checkbox(hs.chrShiftTab,'Provide chromatic shift correction',[x y w h],@chromaticShiftCB,tinyfont);
  y = y-h;
  hs.minChrShiftSpotsText = label(hs.chrShiftTab,'Min spots per jobset',[x y labelw h],tinyfont);
  hs.minChrShiftSpots = editbox(hs.chrShiftTab,[],[editx y editw h],tinyfont);
  
  % Constructing chromatic shift panel
  colwidth = 55;
  x = 2.5;
  w = colwidth-5;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  
  hs.chrShiftPanel = uipanel(hs.chrShiftTab,'Units','characters','Position',[x y-5.75*h w 5.75*h],'FontSize',tinyfont,'Title','Chromatic shift jobsets');
  p = hs.chrShiftPanel;
  y = y-lh;
  labelw = 0.75*w;
  ncols = 4;
  editcol = [(w-2*editw)*(1/ncols), ...
             (w+0.6*editw)*(2/ncols), ...
             (w-0.7*editw)*(3/ncols), ...
             (w-1.5*editw)*(4/ncols)];
  edity = 4.25*h;
  
  % Column headers
  hs.chrShiftChanVectText = label(p,'Channel vector',[1 edity-h labelw/ncols lh],tinyfont);
  hs.chrShiftChanVectText.FontWeight = 'bold';
  hs.chrShiftChanVectText.HorizontalAlignment = 'left';
  hs.chrShiftJobsetText = label(p,'Jobset',[editcol(1) edity-h/2 labelw h],tinyfont);
  hs.chrShiftJobsetText.FontWeight = 'bold';
  hs.chrShiftChanOrderText = label(p,'Channel order',[editcol(2) edity-h labelw/ncols lh],tinyfont);
  hs.chrShiftChanOrderText.FontWeight = 'bold';
  hs.chrShiftChanOrderText.HorizontalAlignment = 'left';
  hs.chrShiftAdjustmentText = label(p,'Adjustment in [x,y,z] (nm)',[editcol(4) edity-h labelw/ncols lh],tinyfont);
  hs.chrShiftAdjustmentText.FontWeight = 'bold';
  hs.chrShiftAdjustmentText.HorizontalAlignment = 'left';
  
  % Jobset, channel order, adjustment, and associated text
  edity = edity-2*h;
  hs.ch1to2Text = label(p,'1  ->  2',[1 edity labelw h],tinyfont);
  hs.ch1to2Arrow = label(p,'->',[(sum(editcol(2:3))+0.75)/2 edity editw h],tinyfont);
  hs.ch1to2 = button(p,'-',[editcol(1) edity 1.9*editw h],@ch1to2CB,tinyfont);
  hs.ch1to2_ch1num = editbox(p,[],[editcol(2) edity editw/5 h],tinyfont);
  hs.ch1to2_ch2num = editbox(p,[],[editcol(3) edity editw/5 h],tinyfont);
  hs.ch1to2_adjText = label(p,'[         ,         ,          ]',[editcol(4) edity editw*2 h],smallfont);
  hs.ch1to2_adjx = editbox(p,[],[editcol(4)+0.75 edity editw/2.7 h],tinyfont);
  hs.ch1to2_adjy = editbox(p,[],[editcol(4)+0.75+editw/2.25 edity editw/2.7 h],tinyfont);
  hs.ch1to2_adjz = editbox(p,[],[editcol(4)+0.75+editw*2/2.25 edity editw/2.7 h],tinyfont);
  edity = edity-h;
  hs.ch1to3Text = label(p,'1  ->  3',[1 edity labelw h],tinyfont);
  hs.ch1to3Arrow = label(p,'->',[(sum(editcol(2:3))+0.75)/2 edity editw h],tinyfont);
  hs.ch1to3 = button(p,'-',[editcol(1) edity 1.9*editw h],@ch1to3CB,tinyfont);
  hs.ch1to3_ch1num = editbox(p,[],[editcol(2) edity editw/5 h],tinyfont);
  hs.ch1to3_ch3num = editbox(p,[],[editcol(3) edity editw/5 h],tinyfont);
  hs.ch1to3_adjText = label(p,'[         ,         ,          ]',[editcol(4) edity editw*2 h],smallfont);
  hs.ch1to3_adjx = editbox(p,[],[editcol(4)+0.75 edity editw/2.7 h],tinyfont);
  hs.ch1to3_adjy = editbox(p,[],[editcol(4)+0.75+editw/2.25 edity editw/2.7 h],tinyfont);
  hs.ch1to3_adjz = editbox(p,[],[editcol(4)+0.75+editw*2/2.25 edity editw/2.7 h],tinyfont);
  edity = edity-h;
  hs.ch2to3Text = label(p,'2  ->  3',[1 edity labelw h],tinyfont);
  hs.ch2to3Arrow = label(p,'->',[(sum(editcol(2:3))+0.75)/2 edity editw h],10);
  hs.ch2to3 = button(p,'-',[editcol(1) edity 1.9*editw h],@ch2to3CB,10);
  hs.ch2to3_ch2num = editbox(p,[],[editcol(2) edity editw/5 h],10);
  hs.ch2to3_ch3num = editbox(p,[],[editcol(3) edity editw/5 h],10);
  hs.ch2to3_adjText = label(p,'[         ,         ,          ]',[editcol(4) edity editw*2 h],smallfont);
  hs.ch2to3_adjx = editbox(p,[],[editcol(4)+0.75 edity editw/2.7 h],tinyfont);
  hs.ch2to3_adjy = editbox(p,[],[editcol(4)+0.75+editw/2.25 edity editw/2.7 h],tinyfont);
  hs.ch2to3_adjz = editbox(p,[],[editcol(4)+0.75+editw*2/2.25 edity editw/2.7 h],tinyfont);
  
  % Chromatic shift filtering options
  colwidth = 35;
  x = 12.5;
  w = colwidth-5;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  y = y-5.5*h;
  hs.chrShiftFilter = checkbox(hs.chrShiftTab,'Filter chromatic shift spots',[x y w h],@chrShiftFilterCB,tinyfont);
  y = y-h;
  hs.chrShiftamplitudeText = label(hs.chrShiftTab,'Min spot intensity (% of max)',[x y labelw h],tinyfont);
  hs.chrShiftamplitude = editbox(hs.chrShiftTab,[],[editx y editw h],10);
  y = y-h;
  hs.chrShiftnnDistText = label(hs.chrShiftTab,'Min spot separation (um)',[x y labelw h],tinyfont);
  hs.chrShiftnnDist = editbox(hs.chrShiftTab,[],[editx y editw h],10);
  y = y-h;
  hs.chrShiftRegion = checkbox(hs.chrShiftTab,'Take only central region',[x y w h],@chrShiftFilterCB,tinyfont);
  y = y-h;
  hs.chrShiftRegionNumText = label(hs.chrShiftTab,'Number of regions (n x n)',[x y labelw h],tinyfont);
  hs.chrShiftRegionNum = editbox(hs.chrShiftTab,[],[editx y editw h],10);
  
  end
  
  %% Spot detection options, tab 2
  
  if tabsToShow(3)
  
  colwidth = 35;
  w = colwidth-5;
  x = 12.5;
  y = toplabely-lh;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  edity = h/4;
  
  hs.gaussFilterSpots = checkbox(hs.spotDetectTab,'Gaussian filter beforehand',[x y w h],[],tinyfont);
  y = y-h;
  t = label(hs.spotDetectTab,'Min spots per frame',[x y labelw h],10);
  hs.minSpotsPerFrame = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  t = label(hs.spotDetectTab,'Max spots per frame',[x y labelw h],10);
  hs.maxSpotsPerFrame = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  % Adaptive spot-detection options
  hs.adaptiveTitle = label(hs.spotDetectTab,'Adaptive',[x y w 1.5],smallfont);
  hs.adaptiveTitle.FontWeight = 'bold';
  y = y-h;
  hs.adaptiveLambdaText = label(hs.spotDetectTab,'Weight for spot count',[x y labelw h],10);
  hs.adaptiveLambda = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  % Neighbour spot-detection options
  hs.neighbourTitle = label(hs.spotDetectTab,'Neighbour',[x y w 1.5],smallfont);
  hs.neighbourTitle.FontWeight = 'bold';
  y = y-h;
  hs.neighbourMaskShapeText = label(hs.spotDetectTab,'Mask shape',[x y labelw h],10);
  hs.neighbourMaskShape = popup(hs.spotDetectTab,neighbourMaskValues,[editx-10 y editw+11 h],@neighbourOptionsCB,10);
  y = y-h;
  hs.neighbourMaskRadiusText = label(hs.spotDetectTab,'Mask radius (um)',[x y labelw h],10);
  hs.neighbourMaskRadius = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  hs.neighbourConeAngleText = label(hs.spotDetectTab,'Cone tip half-angle (deg)',[x y labelw h],tinyfont);
  hs.neighbourConeAngle = editbox(hs.spotDetectTab,[],[editx y editw h],tinyfont);
  hs.neighbourOrientPanel = uipanel(hs.spotDetectTab,'Units','characters','Position',[x y-3*h w 3*h],'FontSize',10,'Title','Neighbour orientation');
  p = hs.neighbourOrientPanel;
  editxl = 10;
  editxr = w-10-editw/3;
  editxc = (editxl+editxr)/2;
  hs.neighbourChanNumText = label(p,'Channel number',[(w-labelw)/2 edity+h labelw h],10);
  hs.neighbourChanNumText.FontWeight = 'bold';
  hs.neighbourChanNumText.HorizontalAlignment = 'center';
  hs.neighbourInnerText = label(p,'inner kchore',[1 edity-h*2/3 labelw/3 2*h],10);
  hs.neighbourOuterText = label(p,'outer kchore',[w-labelw*2/5 edity-h*2/3 labelw/3 2*h],10);
  hs.neighbourOuterText.HorizontalAlignment = 'right';
  hs.neighbourOrient{1} = editbox(p,[],[editxl edity editw/3 h],10);
  hs.neighbourOrient{2} = editbox(p,[],[editxc edity editw/3 h],10);
  hs.neighbourOrient{3} = editbox(p,[],[editxr edity editw/3 h],10);
  % Wavelet spot-detection options
  y = y-4.25*h;
  hs.waveletTitle = label(hs.spotDetectTab,'Wavelet',[x y w 1.5],smallfont);
  hs.waveletTitle.FontWeight = 'bold';
  y = y-h;
  hs.wletLevelThreshText = label(hs.spotDetectTab,'Threshold level for MAD thresholding',[x y-h/2 labelw lh],10);
  hs.wletLevelThresh = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-lh;
  hs.wletLevelAdapt = checkbox(hs.spotDetectTab,'Use adaptive threshold level',[x y w h],[],tinyfont);
  y = y-h;
  hs.wletNumLevelsText = label(hs.spotDetectTab,'Number of levels',[x y labelw h],10);
  hs.wletNumLevels = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  hs.wletLocalMADText = label(hs.spotDetectTab,'Locally-estimated MAD',[x y labelw h],10);
  hs.wletLocalMAD = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  hs.wletBackSubText = label(hs.spotDetectTab,'Background subtraction',[x y labelw h],10);
  hs.wletBackSub = editbox(hs.spotDetectTab,[],[editx y editw h],10);
  y = y-h;
  hs.wletMinLevelText = label(hs.spotDetectTab,'Minimum level',[x y labelw h],10);
  hs.wletMinLevel = editbox(hs.spotDetectTab,[],[editx y editw h],10);

  end
  
  %% Options - MMF options, tab 3
  if tabsToShow(4)
  
  w = colwidth-5;
  x = 12.5;
  y = toplabely-lh;
  lh = 1.5*h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  
  hs.mmfAddSpots = checkbox(hs.mmfTab,'Resolve sub-resolution spots',[x y w h],'',10);
  y = y-h;
  hs.maxMmfTimeText = label(hs.mmfTab,'Max MMF time per frame (s)',[x y labelw h],10);
  hs.maxMmfTime = editbox(hs.mmfTab,[],[editx y editw h],10);
  y = y-h;
  hs.oneBigCluster = checkbox(hs.mmfTab,'Clustering of spots',[x y w h],@mmfClusterCB,10);
  y = y-h;
  hs.clusterSeparationText = label(hs.mmfTab,'Cluster separation (um)',[x y labelw h],10);
  hs.clusterSeparation = editbox(hs.mmfTab,[],[editx y editw h],10);
  y = y-h;
  % form MMF weights panel
  hs.mmfWeightsPanel = uipanel(hs.mmfTab,'Units','characters','Position',[x y-4.5*h w 5.5*h],'FontSize',10,'Title','Restriction weights');
  editw = 0.15*w;
  panx = 1;
  pany = 3.75*h;
  for i=1:3; editx(i) = panx+0.9*w-i*editw-(i-1)*0.25; end
  hs.weightsChannelText{1} = label(hs.mmfWeightsPanel,'Ch.1',[editx(3)+0.5 pany labelw h],10);
  hs.weightsChannelText{2} = label(hs.mmfWeightsPanel,'Ch.2',[editx(2)+0.5 pany labelw h],10);
  hs.weightsChannelText{3} = label(hs.mmfWeightsPanel,'Ch.3',[editx(1)+0.5 pany labelw h],10);
  hs.weightsText.FontWeight = 'bold';
  pany = pany-h;
  hs.alphaAtext = label(hs.mmfWeightsPanel,'Spot intensity',[panx pany labelw h],10);
  for iChan=1:3;
    hs.alphaA{iChan} = editbox(hs.mmfWeightsPanel,[],[editx(4-iChan) pany editw h],10);
  end
  pany = pany-h;
  hs.alphaFtext = label(hs.mmfWeightsPanel,'New spots',[panx pany labelw h],10);
  for iChan=1:3;
    hs.alphaF{iChan} = editbox(hs.mmfWeightsPanel,[],[editx(4-iChan) pany editw h],10);
  end
  pany = pany-h;
  hs.alphaDtext = label(hs.mmfWeightsPanel,'Refined spot location',[panx pany-h/2 labelw/2 lh],10);
  for iChan=1:3;
    hs.alphaD{iChan} = editbox(hs.mmfWeightsPanel,[],[editx(4-iChan) pany editw h],10);
  end
  
  % continue below panel
  y = y-5.75*h;
  editw = 0.2*w;
  editx = x+w-editw;
  hs.mmfTolText = label(hs.mmfTab,'Gaussian fit precision',[x y labelw h],10);
  hs.mmfTol = editbox(hs.mmfTab,[],[editx y editw h],10);
  
  end
  
  %% Options - Tracking options, tab 4
  if tabsToShow(5)
  
  w = colwidth-5;
  lh = 1.5*h;
  y = toplabely-lh;
  
  hs.autoRadii = checkbox(hs.trackingTab,'Calculate search radii from dt',[x y w h],@autoRadiiCB,10);
  y = y-h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  hs.autoRadiidtText = label(hs.trackingTab,'Frame dt',[x y labelw h],10);
  hs.autoRadiidt = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-h;
  hs.autoRadiiAvgDispText = label(hs.trackingTab,'Est. avg. disp. of spots (um/s)',[x y labelw h],10);
  hs.autoRadiiAvgDisp = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-h;
  hs.minSearchRadiusText = label(hs.trackingTab,'Min search radius (um)',[x y labelw h],10);
  hs.minSearchRadius = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-h;
  hs.maxSearchRadiusText = label(hs.trackingTab,'Max search radius (um)',[x y labelw h],10);
  hs.maxSearchRadius = editbox(hs.trackingTab,[],[editx y editw h],10);

  y = y-h;
  hs.useSisterAlignment = checkbox(hs.trackingTab,'Use sister alignment',[x y w h],@useSisterAlignmentCB,10);
  y = y-h;
  % Adjust text box pos for multiple lines.
  hs.maxSisterAlignmentAngleText = label(hs.trackingTab,'Max angle between sisters and plate normal (deg)',[x y-h/2 labelw lh],10);
  hs.maxSisterAlignmentAngle = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-lh;
  hs.maxSisterDistText = label(hs.trackingTab,'Max average distance between sisters (um)',[x y-h/2 labelw lh],10);
  hs.maxSisterDist = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-lh;
  hs.minSisterTrackOverlapText = label(hs.trackingTab,'Min overlap between sister tracks',[x y-h/2 labelw lh],10);
  hs.minSisterTrackOverlap = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-h;
  hs.directionMethodText = label(hs.trackingTab,'Direction assigned by',[x y-h/2 labelw*2/3 lh],10);
  hs.directionMethod = popup(hs.trackingTab,directionMethodValues,[editx-10 y editw+11 h],@directionMethodsCB,tinyfont);
  y = y-h;
  hs.directionWeightText = label(hs.trackingTab,'Weight for direction',[x y-h/2 labelw lh],10);
  hs.directionWeight = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-h;
  hs.directionMinStepsText = label(hs.trackingTab,'Min consecutive time points',[x y-h/2 labelw lh],10);
  hs.directionMinSteps = editbox(hs.trackingTab,[],[editx y editw h],10);
  y = y-h;
  hs.directionSwitchBufferText = label(hs.trackingTab,'Directional switch buffer (time points)',[x y-h/2 labelw lh],10);
  hs.directionSwitchBuffer = editbox(hs.trackingTab,[],[editx y editw h],10);
  
  end
  
  %% Debugging, tab 0
  if tabsToShow(1)
      
      y = toplabely-lh;
      
      hs.asserts = checkbox(hs.debugTab,'Check unexpected problems',[x y labelw h],[],tinyfont);
      y = y-h;
      hs.showPlaneFit = checkbox(hs.debugTab,'Show plane fitting',[x y labelw h],[],tinyfont);
      y = y-h;
      hs.showPlaneFitAll = checkbox(hs.debugTab,'Show for all frames',[x+2.5 y labelw h],[],tinyfont);
%       hs.showPlaneFitText = label(hs.debugTab,'Show plane fitting',[x y labelw h],10);
%       hs.showPlaneFit = popup(hs.debugTab,debugPlaneFitValues,[editx-10 y editw+11 h],[],tinyfont);
      y = y-h;
      hs.gapClosing = checkbox(hs.debugTab,'Show track gap lengths',[x y labelw h],[],tinyfont);
      y = y-h;
      hs.groupSistersText = label(hs.debugTab,'Show sister grouping',[x y labelw h],10);
      hs.groupSisters = popup(hs.debugTab,debugGroupSisterValues,[editx-10 y editw+11 h],[],tinyfont);
      y = y-h;
      hs.showIntensityMasks = checkbox(hs.debugTab,'Show intensity masks',[x y labelw h],[],tinyfont);
      y = y-lh;
      hs.debugMMFtext = label(hs.debugTab,'Mixture model fitting',[x y w 1.5],smallfont);
      hs.debugMMFtext.FontWeight = 'bold';
      y = y-h;
      hs.mmfVerbose = checkbox(hs.debugTab,'Verbose output',[x y labelw h],[],tinyfont);
      y = y-h;
      hs.showMmfClustersText = label(hs.debugTab,'Show clusters',[x y labelw h],10);
      hs.showMmfClusters = popup(hs.debugTab,debugMMFclustersValues,[editx-10 y editw+11 h],[],tinyfont);
      y = y-h;
      hs.showMmfCandsText = label(hs.debugTab,'Show candidates',[x y labelw h],10);
      hs.showMmfCands = popup(hs.debugTab,debugMMFcandsValues,[editx-10 y editw+11 h],[],tinyfont);
      y = y-h;
      hs.showMmfFinalText = label(hs.debugTab,'Show final coordinates',[x y labelw h],10);
      hs.showMmfFinal = popup(hs.debugTab,debugMMFfinalValues,[editx-10 y editw+11 h],[],tinyfont);
      y = y-h;
      hs.showMmfPvals = checkbox(hs.debugTab,'Give p-values',[x y labelw h],[],tinyfont);
      y = y-lh;
      hs.debugCentroidText = label(hs.debugTab,'Centroid fitting',[x y w 1.5],smallfont);
      hs.debugCentroidText.FontWeight = 'bold';
      y = y-h;
      hs.showCentroidFinal = checkbox(hs.debugTab,'Show final coordinates',[x y labelw h],[],tinyfont);
      y = y-lh;
      hs.debugWaveletText = label(hs.debugTab,'Wavelet spot detection',[x y w 1.5],smallfont);
      hs.debugWaveletText.FontWeight = 'bold';
      y = y-h;
      hs.showWavelet = checkbox(hs.debugTab,'Verbose output',[x y labelw h],@showWaveletCB,tinyfont);
      y = y-h;
      hs.saveWaveletImages = checkbox(hs.debugTab,'Save images of output',[x+2.5 y labelw h],[],tinyfont);
      y = y-h;
      hs.showWaveletAdapt = checkbox(hs.debugTab,'Show adaptation of wavelet threshold',[x y 1.5*labelw h],[],tinyfont);
      y = y-lh;
      hs.debugAdaptText = label(hs.debugTab,'Adaptive spot detection',[x y w 1.5],smallfont);
      hs.debugAdaptText.FontWeight = 'bold';
      y = y-h;
      hs.showAdaptive = checkbox(hs.debugTab,'Verbose output',[x y labelw h],[],tinyfont);

  end

  movegui(hs.fig,'center');
  
end

%% Update control status based on contents of jobset.
function updateControls(opts)
  hs = handles;
  
  idx = find(jobset.filename=='/',1,'last');
  if isempty(idx)
    idx=0;
  end
  hs.jobsetName.String = jobset.filename(idx+1:end-4);
  
  % Tracking, tab 4
  hs.autoRadiidt.String = num2str(opts.autoRadiidt);
  hs.autoRadii.Value = ~isempty(opts.autoRadiidt);
  if isempty(opts.autoRadiidt)
    hs.autoRadii.Value = 0; % Off
    hs.autoRadiidt.Enable = 'off';
    hs.autoRadiiAvgDisp.Enable = 'off';
    hs.minSearchRadius.Enable = 'on';
    hs.maxSearchRadius.Enable = 'on';
  else
    hs.autoRadii.Value = 1; % On
    hs.autoRadiidt.Enable = 'on';
    hs.autoRadiiAvgDisp.Enable = 'on';
    hs.minSearchRadius.Enable = 'off';
    hs.maxSearchRadius.Enable = 'off';
  end
  % Assume mean absolute displacment of sisters is about 0.06 μm/s as default.
  hs.autoRadiiAvgDisp.String = num2str(0.06);
  hs.minSearchRadius.String = num2str(opts.minSearchRadius(1));
  hs.maxSearchRadius.String = num2str(opts.maxSearchRadius(1));
  hs.useSisterAlignment.Value = opts.useSisterAlignment;
  hs.maxSisterAlignmentAngle.String = num2str(opts.maxSisterAlignmentAngle);
  if ~hs.useSisterAlignment.Value
    hs.maxSisterAlignmentAngle.Enable = 'off';
  end
  hs.maxSisterDist.String = num2str(opts.maxSisterSeparation);
  hs.minSisterTrackOverlap.String = num2str(opts.minSisterTrackOverlap);
  hs.directionMethod.Value = mapStrings(opts.direction.assignMode,directionMethodValuesJS);
  hs.directionWeight.String = num2str(opts.direction.assignExpWeight);
  hs.directionMinSteps.String = num2str(opts.direction.minConsSteps);
  hs.directionSwitchBuffer.String = num2str(opts.direction.switchBuffer);
  
  % Spot detection, tab 2
  hs.gaussFilterSpots.Value = opts.intensity.gaussFilterSpots;
  hs.minSpotsPerFrame.String = num2str(opts.minSpotsPerFrame);
  hs.maxSpotsPerFrame.String = num2str(opts.maxSpotsPerFrame);
  hs.adaptiveLambda.String = num2str(opts.adaptiveLambda);
  hs.wletLevelThresh.String = num2str(opts.wavelet.levelThresh);
  hs.wletLevelAdapt.Value = opts.wavelet.levelAdapt;
  hs.wletNumLevels.String = num2str(opts.wavelet.numLevels);
  hs.wletLocalMAD.String = num2str(opts.wavelet.localMAD);
  hs.wletBackSub.String = num2str(opts.wavelet.backSub);
  hs.wletMinLevel.String = num2str(opts.wavelet.minLevel);
  hs.neighbourMaskShape.Value = mapStrings(opts.neighbourSpots.maskShape,neighbourMaskValuesJS);
  hs.neighbourMaskRadius.String = num2str(opts.neighbourSpots.maskRadius);
  for iOrient = 1:3
    hs.neighbourOrient{iOrient}.String = num2str(opts.neighbourSpots.channelOrientation(iOrient));
  end
  hs.neighbourConeAngle.String = num2str(opts.neighbourSpots.maskConeAngle);

  
  % Mixture model fitting, tab 3
  hs.mmfAddSpots.Value = opts.mmf.addSpots;
  hs.maxMmfTime.String = num2str(opts.mmf.maxMmfTime);
  hs.clusterSeparation.String = num2str(opts.mmf.clusterSeparation);
  hs.oneBigCluster.Value = ~opts.mmf.oneBigCluster;
  for iChan=1:3
    hs.alphaA{iChan}.String = num2str(opts.mmf.alphaA(iChan));
    hs.alphaD{iChan}.String = num2str(opts.mmf.alphaD(iChan));
    hs.alphaF{iChan}.String = num2str(opts.mmf.alphaF(iChan));
  end
  hs.mmfTol.String = num2str(opts.mmf.mmfTol);
  
  % Chromatic shift, tab 1
  hs.chromaticShift.Value = any(~cellfun('isempty',opts.chrShift.jobset(:)));
  hs.minChrShiftSpots.String = num2str(opts.chrShift.minSpots);
  if ~isempty(opts.chrShift.jobset{1,2})
    hs.ch1to2.String = opts.chrShift.jobset{1,2};
    hs.ch1to2_ch1num.Enable = 'on';
    hs.ch1to2_ch2num.Enable = 'on';
    hs.ch1to2_ch1num.String = num2str(opts.chrShift.chanOrder{1,2}(1));
    hs.ch1to2_ch2num.String = num2str(opts.chrShift.chanOrder{1,2}(2));
    hs.ch1to2_adjx.Enable = 'on';
    hs.ch1to2_adjy.Enable = 'on';
    hs.ch1to2_adjz.Enable = 'on';
    hs.ch1to2_adjx.String = num2str(opts.chrShift.coordinateAdjustments{1,2}(1));
    hs.ch1to2_adjy.String = num2str(opts.chrShift.coordinateAdjustments{1,2}(2));
    hs.ch1to2_adjz.String = num2str(opts.chrShift.coordinateAdjustments{1,2}(3));
  else
    hs.ch1to2.String = '-';
    hs.ch1to2_ch1num.String = '';
    hs.ch1to2_ch2num.String = '';
    hs.ch1to2_ch1num.Enable = 'off';
    hs.ch1to2_ch2num.Enable = 'off';
    hs.ch1to2_adjx.String = '';
    hs.ch1to2_adjy.String = '';
    hs.ch1to2_adjz.String = '';
    hs.ch1to2_adjx.Enable = 'off';
    hs.ch1to2_adjy.Enable = 'off';
    hs.ch1to2_adjz.Enable = 'off';
  end
  if ~isempty(opts.chrShift.jobset{1,3})
    hs.ch1to3.String = opts.chrShift.jobset{1,3};
    hs.ch1to3_ch1num.Enable = 'on';
    hs.ch1to3_ch3num.Enable = 'on';
    hs.ch1to3_ch1num.String = num2str(opts.chrShift.chanOrder{1,3}(1));
    hs.ch1to3_ch3num.String = num2str(opts.chrShift.chanOrder{1,3}(2));
    hs.ch1to3_adjx.Enable = 'on';
    hs.ch1to3_adjy.Enable = 'on';
    hs.ch1to3_adjz.Enable = 'on';
    hs.ch1to3_adjx.String = num2str(opts.chrShift.coordinateAdjustments{1,3}(1));
    hs.ch1to3_adjy.String = num2str(opts.chrShift.coordinateAdjustments{1,3}(2));
    hs.ch1to3_adjz.String = num2str(opts.chrShift.coordinateAdjustments{1,3}(3));
  else
    hs.ch1to3.String = '-';
    hs.ch1to3_ch1num.String = '';
    hs.ch1to3_ch3num.String = '';
    hs.ch1to3_ch1num.Enable = 'off';
    hs.ch1to3_ch3num.Enable = 'off';
    hs.ch1to3_adjx.String = '';
    hs.ch1to3_adjy.String = '';
    hs.ch1to3_adjz.String = '';
    hs.ch1to3_adjx.Enable = 'off';
    hs.ch1to3_adjy.Enable = 'off';
    hs.ch1to3_adjz.Enable = 'off';
  end
  if ~isempty(opts.chrShift.jobset{2,3})
    hs.ch2to3.String = opts.chrShift.jobset{2,3};
    hs.ch2to3_ch2num.Enable = 'on';
    hs.ch2to3_ch3num.Enable = 'on';
    hs.ch2to3_ch2num.String = num2str(opts.chrShift.chanOrder{2,3}(1));
    hs.ch2to3_ch3num.String = num2str(opts.chrShift.chanOrder{2,3}(2));
    hs.ch2to3_adjx.Enable = 'on';
    hs.ch2to3_adjy.Enable = 'on';
    hs.ch2to3_adjz.Enable = 'on';
    hs.ch2to3_adjx.String = num2str(opts.chrShift.coordinateAdjustments{2,3}(1));
    hs.ch2to3_adjy.String = num2str(opts.chrShift.coordinateAdjustments{2,3}(2));
    hs.ch2to3_adjz.String = num2str(opts.chrShift.coordinateAdjustments{2,3}(3));
  else
    hs.ch2to3.String = '-';
    hs.ch2to3_ch2num.String = '';
    hs.ch2to3_ch3num.String = '';
    hs.ch2to3_ch2num.Enable = 'off';
    hs.ch2to3_ch3num.Enable = 'off';
    hs.ch2to3_adjx.String = '';
    hs.ch2to3_adjy.String = '';
    hs.ch2to3_adjz.String = '';
    hs.ch2to3_adjx.Enable = 'off';
    hs.ch2to3_adjy.Enable = 'off';
    hs.ch2to3_adjz.Enable = 'off';
  end
  hs.chrShiftFilter.Value = opts.chrShift.filtering;
  hs.chrShiftamplitude.String = num2str(opts.chrShift.intensityFilter);
  hs.chrShiftnnDist.String = num2str(opts.chrShift.neighbourFilter);
  hs.chrShiftRegionNum.String = num2str(opts.chrShift.regionFilter);
  if opts.chrShift.regionFilter > 1
    hs.chrShiftRegion.Enable = 'on';
    hs.chrShiftRegion.Value = 1;
  else
    hs.chrShiftRegion.Enable = 'off';
    hs.chrShiftRegion.Value = 0;
  end
  
    % general
    hs.asserts.Value = opts.debug.asserts;
    if opts.debug.showPlaneFit == 2
      hs.showPlaneFit.Value = 1;
      hs.showPlaneFitAll.Value = 1;
    else
      hs.showPlaneFit.Value = opts.debug.showPlaneFit;
    end
    hs.groupSisters.Value = opts.debug.groupSisters+1;
    hs.gapClosing.Value = opts.debug.gapClosing;
    hs.showIntensityMasks.Value = opts.debug.showIntensityMasks;
    % centroid
    hs.showCentroidFinal.Value = opts.debug.showCentroidFinal;
    % wavelet
    if opts.debug.showWavelet == 2
      hs.showWavelet.Value = 1;
      hs.saveWaveletImages.Value = 1;
    else
      hs.showWavelet.Value = opts.debug.showWavelet;
    end
    if tabsToShow(1)
      showWaveletCB();
    end
    hs.showWaveletAdapt.Value = opts.debug.showWaveletAdapt;
    % adaptive
    hs.showAdaptive.Value = opts.debug.showAdaptive;
    % MMF
    hs.mmfVerbose.Value = opts.debug.mmfVerbose;
    hs.showMmfCands.Value = find([0 1 -1]==opts.debug.showMmfCands);
    hs.showMmfClusters.Value = find([0 1 -1 -2]==opts.debug.showMmfClusters);
    hs.showMmfFinal.Value = find([0 1 -1]==opts.debug.showMmfFinal);
    hs.showMmfPvals.Value = opts.debug.showMmfPvals;
  
  if ~any(cellfun(@(x) strcmp(x,'neighbour'),opts.spotMode))
    % disable all neighbour options
    hs.neighbourTitle.Enable = 'off';
  end
  if ~any(cellfun(@(x) strcmp(x,'wavelet'),opts.spotMode))
    % disable all wavelet options
    hs.waveletTitle.Enable = 'off';
    hs.wletTitle.Enable = 'off';
    hs.wletLevelThresh.Enable = 'off';
    hs.wletLevelThreshText.Enable = 'off';
    hs.wletLevelAdapt.Enable = 'off';
    hs.wletNumLevels.Enable = 'off';
    hs.wletNumLevelsText.Enable = 'off';
    hs.wletLocalMAD.Enable = 'off';
    hs.wletLocalMADText.Enable = 'off';
    hs.wletBackSub.Enable = 'off';
    hs.wletBackSubText.Enable = 'off';
    hs.wletMinLevel.Enable = 'off';
    hs.wletMinLevelText.Enable = 'off';
    hs.debugWaveletText.Enable = 'off';
    hs.showWavelet.Enable = 'off';
    hs.saveWaveletImages.Enable = 'off';
    hs.showWaveletAdapt.Enable = 'off';
  end
  if ~any(cellfun(@(x) strcmp(x,'adaptive'),opts.spotMode))
    % disable all adaptive options
    hs.adaptiveTitle.Enable = 'off';
    hs.adaptiveLambda.Enable = 'off';
    hs.adaptiveLambdaText.Enable = 'off';
    hs.debugAdaptText.Enable = 'off';
    hs.showAdaptive.Enable = 'off';
  end
  if ~any(cellfun(@(x) strcmp(x,'gaussian'),opts.coordMode))
    % disable all MMF debug options
    hs.debugMMFtext.Enable = 'off';
    hs.mmfVerbose.Enable = 'off';
    hs.showMmfClusters.Enable = 'off';
    hs.showMmfClustersText.Enable = 'off';
    hs.showMmfCands.Enable = 'off';
    hs.showMmfCandsText.Enable = 'off';
    hs.showMmfFinal.Enable = 'off';
    hs.showMmfFinalText.Enable = 'off';
    hs.showMmfPvals.Enable = 'off';
  else
    for iChan=1:3
      if strcmp(opts.coordMode{iChan},'gaussian')
        hs.alphaA{iChan}.Enable = 'on';
        hs.alphaD{iChan}.Enable = 'on';
        hs.alphaF{iChan}.Enable = 'on';
      else
        hs.alphaA{iChan}.Enable = 'off';
        hs.alphaD{iChan}.Enable = 'off';
        hs.alphaF{iChan}.Enable = 'off';
      end
    end
  end
  if ~any(cellfun(@(x) strcmp(x,'centroid'),opts.coordMode))
    % disable centroid debug option
    hs.debugCentroidText.Enable = 'off';
    hs.showCentroidFinal.Enable = 'off';
  end
  if ~strcmp(opts.jobProcess,'zandt')
    % disable all tracking debug options
    hs.gapClosing.Enable = 'off';
    hs.groupSistersText.Enable = 'off';
    hs.groupSisters.Enable = 'off';
  end
  if strcmp(opts.coordSystem,'com')
    % disable plate fitting debug option
    hs.showPlaneFit.Enable = 'off';
    hs.showPlaneFitAll.Enable = 'off';
  end

  handles = hs;
  
  if all(tabsToShow([1,4]))
    mmfClusterCB();
  end
  if tabsToShow(2)
    chromaticShiftCB();
    chrShiftFilterCB();
  end
  if tabsToShow(3)
    neighbourOptionsCB();
  end
  autoRadiiCB();
  if tabsToShow(5)
    directionMethodsCB();
  end
  
end

% Check controls for consistent input.
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

function openExistingCB(hObj,event)
  
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
    r = computeSearchRadii(str2double(handles.autoRadiidt.String),str2double(handles.autoRadiiAvgDisp.String));
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

function chromaticShiftCB(hObj,event)
  if handles.chromaticShift.Value
    handles.minChrShiftSpotsText.Enable = 'on';
    handles.minChrShiftSpots.Enable = 'on';
    handles.chrShiftPanel.ForegroundColor = [0 0 0];
    handles.chrShiftChanVectText.Enable = 'on';
    handles.chrShiftJobsetText.Enable = 'on';
    handles.chrShiftChanOrderText.Enable = 'on';
    handles.chrShiftAdjustmentText.Enable = 'on';
    handles.ch1to2Text.Enable = 'on';
    handles.ch1to2Arrow.Enable = 'on';
    handles.ch1to2.Enable = 'on';
    handles.ch1to3.Enable = 'on';
    handles.ch1to3Text.Enable = 'on';
    handles.ch1to3Arrow.Enable = 'on';
    handles.ch2to3Text.Enable = 'on';
    handles.ch2to3Arrow.Enable = 'on';
    handles.ch2to3.Enable = 'on';
    if all(~strcmp(handles.ch1to2.String,{'-','Unknown source'}))
      handles.ch1to2_ch1num.Enable = 'on';
      handles.ch1to2_ch2num.Enable = 'on';
      handles.ch1to2_adjText.Enable = 'on';
      handles.ch1to2_adjx.Enable = 'on';
      handles.ch1to2_adjy.Enable = 'on';
      handles.ch1to2_adjz.Enable = 'on';
      handles.ch1to2_adjText.Enable = 'on';
    else
      handles.ch1to2_ch1num.String = '';
      handles.ch1to2_ch2num.String = '';
      handles.ch1to2_ch1num.Enable = 'off';
      handles.ch1to2_ch2num.Enable = 'off';
      handles.ch1to2_adjx.Enable = 'off';
      handles.ch1to2_adjy.Enable = 'off';
      handles.ch1to2_adjz.Enable = 'off';
      handles.ch1to2_adjText.Enable = 'off';
    end
    if all(~strcmp(handles.ch1to3.String,{'-','Unknown source'}))
      handles.ch1to3_ch1num.Enable = 'on';
      handles.ch1to3_ch3num.Enable = 'on';
      handles.ch1to3_adjx.Enable = 'on';
      handles.ch1to3_adjy.Enable = 'on';
      handles.ch1to3_adjz.Enable = 'on';
      handles.ch1to3_adjText.Enable = 'on';
    else
      handles.ch1to3_ch1num.String = '';
      handles.ch1to3_ch3num.String = '';
      handles.ch1to3_ch1num.Enable = 'off';
      handles.ch1to3_ch3num.Enable = 'off';
      handles.ch1to3_adjx.Enable = 'off';
      handles.ch1to3_adjy.Enable = 'off';
      handles.ch1to3_adjz.Enable = 'off';
      handles.ch1to3_adjText.Enable = 'off';
    end
    if all(~strcmp(handles.ch2to3.String,{'-','Unknown source'}))
      handles.ch2to3_ch2num.Enable = 'on';
      handles.ch2to3_ch3num.Enable = 'on';
      handles.ch2to3_adjx.Enable = 'on';
      handles.ch2to3_adjy.Enable = 'on';
      handles.ch2to3_adjz.Enable = 'on';
      handles.ch2to3_adjText.Enable = 'on';
    else
      handles.ch2to3_ch2num.String = '';
      handles.ch2to3_ch3num.String = '';
      handles.ch2to3_ch2num.Enable = 'off';
      handles.ch2to3_ch3num.Enable = 'off';
      handles.ch2to3_adjx.Enable = 'off';
      handles.ch2to3_adjy.Enable = 'off';
      handles.ch2to3_adjz.Enable = 'off';
      handles.ch2to3_adjText.Enable = 'off';
    end
    handles.chrShiftFilter.Enable = 'on';
  else
    handles.minChrShiftSpotsText.Enable = 'off';
    handles.minChrShiftSpots.Enable = 'off';
    handles.chrShiftPanel.ForegroundColor = [0.5 0.5 0.5];
    handles.chrShiftChanVectText.Enable = 'off';
    handles.chrShiftJobsetText.Enable = 'off';
    handles.chrShiftChanOrderText.Enable = 'off';
    handles.chrShiftAdjustmentText.Enable = 'off';
    handles.ch1to2Text.Enable = 'off';
    handles.ch1to2Arrow.Enable = 'off';
    handles.ch1to2.Enable = 'off';
    handles.ch1to3Text.Enable = 'off';
    handles.ch1to3Arrow.Enable = 'off';
    handles.ch1to3.Enable = 'off';
    handles.ch2to3Text.Enable = 'off';
    handles.ch2to3Arrow.Enable = 'off';
    handles.ch2to3.Enable = 'off';
    handles.ch1to2_ch1num.Enable = 'off';
    handles.ch1to2_ch2num.Enable = 'off';
    handles.ch1to3_ch1num.Enable = 'off';
    handles.ch1to3_ch3num.Enable = 'off';
    handles.ch2to3_ch2num.Enable = 'off';
    handles.ch2to3_ch3num.Enable = 'off';
    handles.ch1to2_adjText.Enable = 'off';
    handles.ch1to2_adjx.Enable = 'off';
    handles.ch1to2_adjy.Enable = 'off';
    handles.ch1to2_adjz.Enable = 'off';
    handles.ch1to3_adjText.Enable = 'off';
    handles.ch1to3_adjx.Enable = 'off';
    handles.ch1to3_adjy.Enable = 'off';
    handles.ch1to3_adjz.Enable = 'off';
    handles.ch2to3_adjText.Enable = 'off';
    handles.ch2to3_adjx.Enable = 'off';
    handles.ch2to3_adjy.Enable = 'off';
    handles.ch2to3_adjz.Enable = 'off';
    handles.chrShiftFilter.Value = 0;
    chrShiftFilterCB();
    handles.chrShiftFilter.Enable = 'off';
  end
end

function chrShiftFilterCB(hObj,event)
  if handles.chrShiftFilter.Value
    handles.chrShiftamplitudeText.Enable = 'on';
    handles.chrShiftamplitude.Enable = 'on';
    handles.chrShiftnnDistText.Enable = 'on';
    handles.chrShiftnnDist.Enable = 'on';
    handles.chrShiftRegion.Enable = 'on';
  else
    handles.chrShiftamplitudeText.Enable = 'off';
    handles.chrShiftamplitude.Enable = 'off';
    handles.chrShiftnnDistText.Enable = 'off';
    handles.chrShiftnnDist.Enable = 'off';
    handles.chrShiftRegion.Value = 0;
    handles.chrShiftRegion.Enable = 'off';
  end
  if handles.chrShiftRegion.Value
    handles.chrShiftRegionNumText.Enable = 'on';
    handles.chrShiftRegionNum.Enable = 'on';
  else
    handles.chrShiftRegionNumText.Enable = 'off';
    handles.chrShiftRegionNum.String = 1;
    handles.chrShiftRegionNum.Enable = 'off';
  end
end

function ch1to2CB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  handles.ch1to2.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{1,2} = file;
  jobset.options.chrShift.jobset{2,1} = file;
  handles.ch1to2_ch1num.Enable = 'on';
  handles.ch1to2_ch1num.String = jobset.options.chrShift.chanOrder{1,2}(1);
  handles.ch1to2_ch2num.Enable = 'on';
  handles.ch1to2_ch2num.String = jobset.options.chrShift.chanOrder{1,2}(2);
end

function ch1to3CB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  handles.ch1to3.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{1,3} = file;
  jobset.options.chrShift.jobset{3,1} = file;
  handles.ch1to3_ch1num.Enable = 'on';
  handles.ch1to3_ch1num.String = jobset.options.chrShift.chanOrder{1,3}(1);
  handles.ch1to3_ch3num.Enable = 'on';
  handles.ch1to3_ch3num.String = jobset.options.chrShift.chanOrder{1,3}(2);
end

function ch2to3CB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  handles.ch2to3.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{2,3} = file;
  jobset.options.chrShift.jobset{3,2} = file;
  handles.ch2to3_ch2num.Enable = 'on';
  handles.ch2to3_ch2num.String = jobset.options.chrShift.chanOrder{2,3}(1);
  handles.ch2to3_ch3num.Enable = 'on';
  handles.ch2to3_ch3num.String = jobset.options.chrShift.chanOrder{2,3}(2);
end

function neighbourOptionsCB(hObj,event)
  neighChans = cellfun(@(x) strcmp(x,'neighbour'),jobset.options.spotMode);
  switch sum(neighChans)
    case 0
      handles.neighbourMaskShapeText.Enable = 'off';
      handles.neighbourMaskShape.Enable = 'off';
      handles.neighbourMaskRadiusText.Enable = 'off';
      handles.neighbourMaskRadius.Enable = 'off';
      handles.neighbourOrientPanel.ForegroundColor = [0.5 0.5 0.5];
      handles.neighbourChanVectText.Enable = 'off';
      handles.neighbourChanNumText.Enable = 'off';
      handles.neighbourInnerText.Enable = 'off';
      for iChan=1:3
        handles.neighbourOrient{iChan}.Enable = 'off';
      end
      handles.neighbourOuterText.Enable = 'off';
    case 1
      handles.neighbourMaskShapeText.Enable = 'on';
      handles.neighbourMaskShape.Enable = 'on';
      handles.neighbourMaskRadiusText.Enable = 'on';
      handles.neighbourMaskRadius.Enable = 'on';
      handles.neighbourOrientPanel.ForegroundColor = [0 0 0];
      handles.neighbourChanVectText.Enable = 'on';
      handles.neighbourChanNumText.Enable = 'on';
      handles.neighbourInnerText.Enable = 'on';
      handles.neighbourOuterText.Enable = 'on';
      for iChan=1:3
        if ismember(iChan,[jobset.options.coordSystemChannel find(neighChans)])
          handles.neighbourOrient{iChan}.Enable = 'on';
        else
          handles.neighbourOrient{iChan}.Enable = 'off';
        end
      end
    case 2
      handles.neighbourMaskShapeText.Enable = 'on';
      handles.neighbourMaskShape.Enable = 'on';
      handles.neighbourMaskRadiusText.Enable = 'on';
      handles.neighbourMaskRadius.Enable = 'on';
      handles.neighbourOrientPanel.ForegroundColor = [0 0 0];
      handles.neighbourChanVectText.Enable = 'on';
      handles.neighbourChanNumText.Enable = 'on';
      handles.neighbourInnerText.Enable = 'on';
      handles.neighbourOuterText.Enable = 'on';
      for iChan = 1:3
        handles.neighbourOrient{iChan}.Enable = 'on';
      end
  end

  if strcmp(mapStrings(handles.neighbourMaskShape.Value,neighbourMaskValues),'Circle')
    handles.neighbourOrientPanel.ForegroundColor = [0.5 0.5 0.5];
      handles.neighbourChanVectText.Enable = 'off';
      handles.neighbourChanNumText.Enable = 'off';
      handles.neighbourInnerText.Enable = 'off';
      handles.neighbourOuterText.Enable = 'off';
      for iChan=1:3
        handles.neighbourOrient{iChan}.Enable = 'off';
      end
  end    
  if strcmp(mapStrings(handles.neighbourMaskShape.Value,neighbourMaskValues),'Cone')
    handles.neighbourConeAngleText.Enable = 'on';
    handles.neighbourConeAngle.Enable = 'on';
  else
    handles.neighbourConeAngleText.Enable = 'off';
    handles.neighbourConeAngle.Enable = 'off';
  end
  
  if strcmp(jobset.options.jobProcess,'chrShift') || ...
          strcmp(jobset.options.coordSystem,'com')
    handles.neighbourMaskShapeText.Enable = 'off';
    handles.neighbourMaskShape.Value = 1;
    handles.neighbourMaskShape.Enable = 'off';
    handles.neighbourOrientPanel.ForegroundColor = [0.5 0.5 0.5];
    handles.neighbourChanVectText.Enable = 'off';
    handles.neighbourChanNumText.Enable = 'off';
    handles.neighbourInnerText.Enable = 'off';
    handles.neighbourOuterText.Enable = 'off';
    for iChan=1:3
        handles.neighbourOrient{iChan}.Enable = 'off';
    end
  end
end

function mmfClusterCB(hObj,event)
  if handles.oneBigCluster.Value
    handles.clusterSeparationText.Enable = 'on';
    handles.clusterSeparation.Enable = 'on';
  else
    handles.clusterSeparationText.Enable = 'off';
    handles.clusterSeparation.Enable = 'off';
  end
end

function directionMethodsCB(hObj,event)
  if strcmp(mapStrings(handles.directionMethod.Value,directionMethodValues),'Voting')
    handles.directionWeight.Enable = 'on';
    handles.directionWeightText.Enable = 'on';
    handles.directionMinSteps.Enable = 'off';
    handles.directionMinStepsText.Enable = 'off';
    handles.directionSwitchBuffer.Enable = 'off';
    handles.directionSwitchBufferText.Enable = 'off';
  elseif strcmp(mapStrings(handles.directionMethod.Value,directionMethodValues),'Step-wise')
    handles.directionWeight.Enable = 'off';
    handles.directionWeightText.Enable = 'off';
    handles.directionMinSteps.Enable = 'on';
    handles.directionMinStepsText.Enable = 'on';
    handles.directionSwitchBuffer.Enable = 'on';
    handles.directionSwitchBufferText.Enable = 'on';
  end
end

function showWaveletCB(hObj,event)
  if handles.showWavelet.Value
    handles.saveWaveletImages.Enable = 'on';
  else
    handles.saveWaveletImages.Enable = 'off';
    handles.saveWaveletImages.Value = 0;
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

function newJobsetCB(hObj,event)
  newJobsetName = [jobset.movieDirectory '/' handles.jobsetName.String '.mat'];
  if ~strcmp(jobset.filename,newJobsetName)
    if ~checkControls()
      return
    end
    updateJobset();
    jobset.filename = newJobsetName;
    kitSaveJobset(jobset);
    uiresume(gcf);
  else
    errorbox('Jobset name is currently unchanged. Change the name to save a new jobset, otherwise instead click ''Save''.');
    return
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

  opts = jobset.options;
  
  % tracking
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
  direction = opts.direction;
  direction.assignMode = mapStrings(handles.directionMethod.Value,directionMethodValuesJS);
  direction.assignExpWeight = str2double(handles.directionWeight.String);
  direction.minConsSteps = str2double(handles.directionMinSteps.String);
  direction.switchBuffer = str2double(handles.directionSwitchBuffer.String);
  opts.direction = direction;
  
  % spot detection
  opts.intensity.gaussFilterSpots = handles.gaussFilterSpots.Value;
  opts.minSpotsPerFrame = str2double(handles.minSpotsPerFrame.String);
  opts.maxSpotsPerFrame = str2double(handles.maxSpotsPerFrame.String);
  opts.adaptiveLambda = str2double(handles.adaptiveLambda.String);
  wavelet = opts.wavelet;
  wavelet.levelThresh = str2double(handles.wletLevelThresh.String);
  wavelet.levelAdapt = handles.wletLevelAdapt.Value;
  wavelet.numLevels = str2double(handles.wletNumLevels.String);
  wavelet.localMAD = str2double(handles.wletLocalMAD.String);
  wavelet.backSub = str2double(handles.wletBackSub.String);
  wavelet.minLevel = str2double(handles.wletMinLevel.String);
  opts.wavelet = wavelet;
  neighbourSpots = opts.neighbourSpots;
  neighbourSpots.maskShape = mapStrings(handles.neighbourMaskShape.Value,neighbourMaskValuesJS);
  neighbourSpots.maskRadius = str2double(handles.neighbourMaskRadius.String);
  neighbourSpots.maskConeAngle = str2double(handles.neighbourConeAngle.String);
  for iOrient=1:3;
    neighbourSpots.channelOrientation(iOrient) = str2double(handles.neighbourOrient{iOrient}.String);
  end
  opts.neighbourSpots = neighbourSpots;
  
  % mmf
  mmf  = opts.mmf;
  mmf.addSpots = handles.mmfAddSpots.Value;
  mmf.maxMmfTime = str2double(handles.maxMmfTime.String);
  mmf.clusterSeparation = str2double(handles.clusterSeparation.String);
  mmf.oneBigCluster = ~handles.oneBigCluster.Value;
  for iChan = 1:3;
    mmf.alphaA(iChan) = str2double(handles.alphaA{iChan}.String);
    mmf.alphaD(iChan) = str2double(handles.alphaD{iChan}.String);
    mmf.alphaF(iChan) = str2double(handles.alphaF{iChan}.String);
  end
  mmf.mmfTol = str2double(handles.mmfTol.String);
  opts.mmf = mmf;
  
  % chromatic shift
  if handles.chromaticShift.Value
    chrShift = opts.chrShift;
    chrShift.minSpots = str2double(handles.minChrShiftSpots.String);
    % channel orders
    chrShift.chanOrder{1,2}(1) = str2double(handles.ch1to2_ch1num.String);
    chrShift.chanOrder{1,2}(2) = str2double(handles.ch1to2_ch2num.String);
    chrShift.chanOrder{2,1} = fliplr(chrShift.chanOrder{1,2});
    chrShift.chanOrder{1,3}(1) = str2double(handles.ch1to3_ch1num.String);
    chrShift.chanOrder{1,3}(2) = str2double(handles.ch1to3_ch3num.String);
    chrShift.chanOrder{3,1} = fliplr(chrShift.chanOrder{1,3});
    chrShift.chanOrder{2,3}(1) = str2double(handles.ch2to3_ch2num.String);
    chrShift.chanOrder{2,3}(2) = str2double(handles.ch2to3_ch3num.String);
    chrShift.chanOrder{3,2} = fliplr(chrShift.chanOrder{2,3});
    % any adjustments
    chrShift.coordinateAdjustments{1,2}(1) = str2double(handles.ch1to2_adjx.String);
    chrShift.coordinateAdjustments{1,2}(2) = str2double(handles.ch1to2_adjy.String);
    chrShift.coordinateAdjustments{1,2}(3) = str2double(handles.ch1to2_adjz.String);
    chrShift.coordinateAdjustments{2,1} = -chrShift.coordinateAdjustments{1,2};
    chrShift.coordinateAdjustments{1,3}(1) = str2double(handles.ch1to3_adjx.String);
    chrShift.coordinateAdjustments{1,3}(2) = str2double(handles.ch1to3_adjy.String);
    chrShift.coordinateAdjustments{1,3}(3) = str2double(handles.ch1to3_adjz.String);
    chrShift.coordinateAdjustments{3,1} = -chrShift.coordinateAdjustments{1,3};
    chrShift.coordinateAdjustments{2,3}(1) = str2double(handles.ch2to3_adjx.String);
    chrShift.coordinateAdjustments{2,3}(2) = str2double(handles.ch2to3_adjy.String);
    chrShift.coordinateAdjustments{2,3}(3) = str2double(handles.ch2to3_adjz.String);
    chrShift.coordinateAdjustments{3,2} = -chrShift.coordinateAdjustments{2,3};
    % filtering
    opts.chrShift.filtering = handles.chrShiftFilter.Value;
    if handles.chrShiftFilter.Value
      chrShift.filtering = 1;
      chrShift.intensityFilter = str2double(handles.chrShiftamplitude.String);
      chrShift.neighbourFilter = str2double(handles.chrShiftnnDist.String);
      if handles.chrShiftRegion.Value
        chrShift.regionFilter = str2double(handles.chrShiftRegionNum.String);
      else
        chrShift.regionFilter = 1;
      end
    else
      chrShift.intensityFilter = 0;
      chrShift.neighbourFilter = 0;
      chrShift.regionFilter = 1;
    end
    result = getChromaticShiftResults(chrShift);
    chrShift.result = result;
    opts.chrShift = chrShift;
  end
  
  if tabsToShow(1)
    % general
    opts.debug.asserts = handles.asserts.Value;
    opts.debug.showPlaneFit = handles.showPlaneFit.Value + handles.showPlaneFitAll.Value;
    opts.debug.groupSisters = handles.groupSisters.Value-1;
    opts.debug.gapClosing = handles.gapClosing.Value;
    opts.debug.showIntensityMasks = handles.showIntensityMasks.Value;
    % centroid
    opts.debug.showCentroidFinal = handles.showCentroidFinal.Value;
    % wavelet
    opts.debug.showWavelet = handles.showWavelet.Value + handles.saveWaveletImages.Value;
    opts.debug.showWaveletAdapt = handles.showWaveletAdapt.Value;
    % adaptive
    opts.debug.showAdaptive = handles.showAdaptive.Value;
    % MMF
    opts.debug.mmfVerbose = handles.mmfVerbose.Value;
    debugCodes = [0 1 -1];
    opts.debug.showMmfCands = debugCodes(handles.showMmfCands.Value);
    opts.debug.showMmfClusters = debugCodes(handles.showMmfClusters.Value);
    opts.debug.showMmfFinal = debugCodes(handles.showMmfFinal.Value);
    opts.debug.showMmfPvals = handles.showMmfPvals.Value;
  end
  
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
	  if i>=j || any(strcmp(jS,{'Unknown source','-'})); continue; end
      jS = kitLoadJobset(jS);
      neighFilt = chrShift.neighbourFilter;
      intFilt = chrShift.intensityFilter/100;
      regFilt = chrShift.regionFilter;
      mS = kitLoadAllJobs(jS);
      if chrShift.filtering
        for iMov = 1:length(mS)
          mS{iMov} = chrsFilterSpots(mS{iMov}, ...
              'neighbourFilter',neighFilt,'intensityFilter',intFilt,'regionFilter',regFilt, ...
              'revert',1,'referenceChan',handles.coordSysChNum);
        end
      end
      [result,~] = chrsCalculateChromaticShift(mS,[i j],...
          'filtered',chrShift.filtering,'interphaseToMetaphase',chrShift.interphase);
      cellResult{i,j} = result; cellResult{j,i} = result.*[-1 -1 -1 1 1 1];
    end
  end       
end

end % kitUpdateOptions
