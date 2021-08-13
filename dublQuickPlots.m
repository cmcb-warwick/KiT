function dublQuickPlots()
% DUBLQUICKPLOTS Displays GUI to allow user to produce plots of various
% dual-channel analyses.
%
% When loading experiments, need to load saved output from
% dublIntraMeasurements.
%
% Copyright (c) 2017 C. A. Smith

% Set up GUI.
dockStatus = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal');
handles = createControls();
updateControls();
handles.fig.Visible = 'on';
uiwait(gcf);
set(0,'DefaultFigureWindowStyle',dockStatus);

%% NESTED FUNCTIONS

% Create all main controls.
function hs = createControls()
    
  % give an empty dataset for now
  hs.data = {[]};
    
  % Create figure.
  figw = 90;
  figh = 29;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = 'Dual-channel plotting tools';
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';
  
  % Define font sizes.
  medfont = 14;
  smallfont = 12;
  tinyfont = 10;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  toplabely = figh; %top-most point
  
  % Experiments.
  y = toplabely-lh; x = dx; txtw = 20;
  t = label(hs.fig,'Experiments',[x y txtw h],medfont);
  t.FontWeight = 'bold';
  
  % Add experiment.
  btnw = 5; x = x+txtw+dx;
  hs.addexptBtn = button(hs.fig,'Add',[x y btnw h],@addexptCB);
  % List structure for added experiments.
  btnw = 2; labw = 6; txtw = 30;
  y = y-lh;
  for iExpt = 1:5
    x = figw/2-(txtw+labw+2*btnw+3*dx);
    hs.exptFile{iExpt} = label(hs.fig,'',[x y txtw h],smallfont);
    hs.exptFile{iExpt}.HorizontalAlignment = 'right';
    x = figw/2-(labw+2*btnw+2*dx);
    hs.exptLabel{iExpt} = editbox(hs.fig,[],[x y labw h],tinyfont);
    hs.exptLabel{iExpt}.Enable = 'off';
    x = figw/2-(2*btnw+dx);
    hs.pairedTxt{iExpt} = label(hs.fig,'',[x y labw h],tinyfont);
    hs.pairedTxt{iExpt}.Enable = 'off';
    x = figw/2-(btnw+dx);
    hs.rmvexptBtn{iExpt} = button(hs.fig,'-',[x y btnw h*3/4],@rmvexptCB,tinyfont);
    hs.rmvexptBtn{iExpt}.Enable = 'off';
    y = y-h;
  end
  hs.nExpts = 0;
  y = y+h;
  
  % Change reference frame here.
  y = y-lh; x = dx; ddx=0.5; radw = 8;
  txtw = figw/2-(3*radw+2*ddx+3*dx);
  refFrame = {'sisters','plate','image'};
  label(hs.fig,'Reference frame',[x y txtw lh]);
  x = x+(txtw+dx);
  for i=1:3
    hs.refFrame{i} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String',refFrame{i},'Position',[x y radw h],'Callback',@refFrameCB);
    x = x+(radw+ddx);
  end
  hs.refFrame{3}.Value = 1;
  hs.refFrameVal = 'microscope';
  
  % Filtering option.
  y = y-h;
  txtw = 28; x = dx; 
  hs.filtering = checkbox(hs.fig,'Apply filter in z-directional delta',[x y txtw h],'',smallfont);
  hs.filtering.Value = 1;
  
  % Button widths and separations for all measurements.
  btnw = 3.5; ddx = 1;
  
  % Measurements.
  y = y-lh; x = dx; txtw = 20;
  t = label(hs.fig,'Measurements',[x y txtw h],medfont);
  t.FontWeight = 'bold';

  % K-K.
  y = y-h;
  x = dx; nBtns = 5; meas={'<html>d<sub>x</html>',...
                           '<html>d<sub>y</html>',...
                           '<html>d<sub>z</html>',...
                           '<html>d<sub>2D</html>',...
                           '<html>d<sub>3D</html>'};
  txtw = figw/2-(3*dx + nBtns*btnw + (nBtns-1)*ddx);
  hs.kkText = label(hs.fig,'Sister separation:',[x y txtw h],smallfont);
  x = x+(txtw+dx);
  for iMeas = 1:nBtns
    hs.kkBtn{iMeas} = button(hs.fig,meas{iMeas},[x y btnw h],@plotKKCB,tinyfont);
    x = x+(btnw+ddx);
  end
  
  % Delta.
  y = y-lh;
  x = dx; nBtns = 6; meas={'<html>&Delta<sub>x</html>',...
                           '<html>&Delta<sub>y</html>',...
                           '<html>&Delta<sub>z</html>',...
                           '<html>&Delta<sub>1D</html>',...
                           '<html>&Delta<sub>2D</html>',...
                           '<html>&Delta<sub>3D</html>'};
  txtw = figw/2-(3*dx + nBtns*btnw + (nBtns-1)*ddx);
  hs.delText = label(hs.fig,'Delta:',[x y txtw h],smallfont);
  x = x+(txtw+dx);
  for iMeas = 1:nBtns
    hs.delBtn{iMeas} = button(hs.fig,meas{iMeas},[x y btnw h],@plotDelCB,tinyfont);
    x = x+(btnw+ddx);
  end
  
  % Swivel.
  y = y-lh;
  x = dx; nBtns = 3; meas={'<html>&theta<sub>y</html>',...
                           '<html>&theta<sub>z</html>',...
                           '<html>&theta<sub>3D</html>'};
  txtw = figw/2-(3*dx + nBtns*btnw + (nBtns-1)*ddx);
  hs.swivText = label(hs.fig,'Swivel:',[x y txtw h],smallfont);
  x = x+(txtw+dx);
  for iMeas = 1:nBtns
    hs.swivBtn{iMeas} = button(hs.fig,meas{iMeas},[x y btnw h],@plotSwivCB,tinyfont);
    x = x+(btnw+ddx);
  end
  
  % Correlations.
  y = y-lh;
  btnw = 7; x = dx; txtw = 20;
  hs.corrText = label(hs.fig,'Correlations',[x y txtw h],medfont);
  hs.corrText.FontWeight = 'bold';
  
  % 2D and 3D.
  y = y-h;
  nBtns = 2; meas={'<html>&Delta vs d</html>',...
                   '<html>&Delta vs &theta</html>'};
  txtw = figw/2-(3*dx + nBtns*btnw + (nBtns-1)*ddx);
  hs.corr2dText = label(hs.fig,'2D:',[x y 5 h],smallfont);
  x = x+(txtw+dx);
  for iMeas = 1:nBtns
    hs.corr2dBtn{iMeas} = button(hs.fig,meas{iMeas},[x y btnw h],@plot2DcorrCB,tinyfont);
    x = x+(btnw+ddx);
  end
  y = y-h; x = dx;
  hs.corr3dText = label(hs.fig,'3D:',[x y 5 h],smallfont);
  x = x+(txtw+dx);
  for iMeas = 1:nBtns
    hs.corr3dBtn{iMeas} = button(hs.fig,meas{iMeas},[x y btnw h],@plot3DcorrCB,tinyfont);
    x = x+(btnw+ddx);
  end
  
  % Create figure panel.
  sfigw = (figw-3*dx)/2; sfigh = 20;
  sfigy = toplabely-(sfigh+dx/4);
  x = figw/2+dx;
  hs.figPanel = uipanel(hs.fig,'Units','characters',...
      'Position',[x sfigy sfigw sfigh],'FontSize',smallfont,'BackgroundColor','w');
  hs.subfig = subplot(1,1,1,'Parent',hs.figPanel);
  
  % Figure export and save buttons.
  y = sfigy-(2*h); x = figw/2+dx;
  ddx = 1;
  btnw = 12;
  hs.savebtn = button(hs.fig,'Save figure',[x y btnw lh],@saveCB,smallfont);
  hs.savebtn.Enable = 'off'; % silence until plot produced
  x = x+btnw+ddx;
  hs.exportbtn = button(hs.fig,'Export figure',[x y btnw lh],@exportCB,smallfont);
  hs.exportbtn.Enable = 'off'; % silence until plot produced
  
  % 'Close' button.
  btnw = 10;
  x = figw-(btnw+dx);
  y = 1;
  hs.closebtn = button(hs.fig,'Close',[x y btnw lh],@closeCB,medfont);
  
  movegui(hs.fig,'center');
  
end

function updateControls()
     
    hs = handles;
    
    % check whether have any experiments without paired data
    if hs.nExpts == 0
        paired = 0;
    else
        % loop over experiments - assume paired until proven otherwise
        paired=1; iExpt=0;
        while iExpt<hs.nExpts && paired
            iExpt=iExpt+1;
            paired = strcmp(hs.pairedTxt{iExpt}.String,'p');
        end
    end
    hs.paired = paired;
    
    % update environment as required
    if paired
%         hs.refFrame{1}.Enable = 'on'; % SUPPRESS THIS FOR NOW
        for i=1:5
            hs.kkBtn{i}.Enable = 'on';
        end
        hs.delBtn{4}.Enable = 'on';
        for i=1:3
            hs.swivBtn{i}.Enable = 'on';
        end
        for i=1:2
            hs.corr2dBtn{i}.Enable = 'on';
            hs.corr3dBtn{i}.Enable = 'on';
        end
    else
        hs.refFrame{1}.Enable = 'off';
        for i=1:5
            hs.kkBtn{i}.Enable = 'off';
        end
        hs.delBtn{4}.Enable = 'off';
        for i=1:3
            hs.swivBtn{i}.Enable = 'off';
        end
        for i=1:2
            hs.corr2dBtn{i}.Enable = 'off';
            hs.corr3dBtn{i}.Enable = 'off';
        end
    end
    
    % can only plot correlations for a single experiment
    if hs.nExpts > 1
        for i=1:2
            hs.corr2dBtn{i}.Enable = 'off';
            hs.corr3dBtn{i}.Enable = 'off';
        end
    end
    
    handles = hs;
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

function closeCB(hObj,event)
  % close the figure 
  close(handles.fig);
end

function addexptCB(hObj,event)
  % ask user to get an intraMeasurements file
  [filename,pathname] = uigetfile('*.mat','Select an intra measurements structure file');
  filepath = fullfile(pathname,filename);
  
  % use the filepath to get the data
  nExpts = handles.nExpts+1;
  data = loadData(filepath);
  handles.data{nExpts} = data;
  
  % change text to the file name
  maxDirLen = 28;
  filename = filename(1:end-4);
  filename = strshorten(filename,maxDirLen);
  set(handles.exptFile{nExpts},'String',filename);
  
  % check whether the data is paired
  handles.pairedTxt{nExpts}.Enable = 'on';
  if isfield(data.microscope,'sisSep') && ~isempty(data.microscope.sisSep.threeD)
    set(handles.pairedTxt{nExpts},'String','p')
    set(handles.pairedTxt{nExpts},'TooltipString', 'This file contains paired data');
  else
    set(handles.pairedTxt{nExpts},'String','')
    set(handles.pairedTxt{nExpts},'TooltipString', '');
  end
  
  % change environment accordingly
  handles.exptLabel{nExpts}.Enable = 'on';
  handles.rmvexptBtn{nExpts}.Enable = 'on';
  
  % reached capacity
  if nExpts == 5
      handles.addexptBtn.Enable = 'off';
  end
  handles.nExpts = nExpts;
  
  % update based on pairedness
  updateControls;
  
end

function rmvexptCB(hObj,event)
  
  % get number of expts
  nExpts = handles.nExpts;
  
  % find which experiment needs removing
  ypos = hObj.Position(2);
  remidx = 0; stopd = 0;
  while ~stopd
      remidx = remidx+1;
      stopd = (handles.rmvexptBtn{remidx}.Position(2) == ypos);
  end
  
  % change text to the file name
  handles.data(remidx) = [];
  for iExpt = remidx:nExpts-1
    filename = handles.exptFile{iExpt+1}.String;
    set(handles.exptFile{iExpt},'String',filename);
    exptlabel = handles.exptLabel{iExpt+1}.String;
    set(handles.exptLabel{iExpt},'String',exptlabel);
    pairlabel = handles.pairedTxt{iExpt+1}.String;
    set(handles.pairedTxt{iExpt},'String',pairlabel);
    pairttstring = handles.pairedTxt{iExpt+1}.TooltipString;
    set(handles.pairedTxt{iExpt},'TooltipString', pairttstring);
  end
  
  % change environment accordingly
  set(handles.exptFile{nExpts},'String','');
  set(handles.exptLabel{nExpts},'String','');
  handles.exptLabel{nExpts}.Enable = 'off';
  handles.rmvexptBtn{nExpts}.Enable = 'off';
  handles.pairedTxt{nExpts}.String = '';
  handles.pairedTxt{nExpts}.TooltipString = '';
  
  % reached capacity
  handles.nExpts = nExpts-1;
  
  % update paired status
  updateControls;
  
end

function refFrameCB(hObj,event)
    
    if exist('hObj')
        rfoptions = {'sisters','plate','image'};
        % get which radio button was selected
        refFrameVal = hObj.String;
        idx = find(strcmp(refFrameVal,rfoptions));
        % untick all other buttons
        for notIdx = setdiff(1:3,idx)
          handles.refFrame{notIdx}.Value = 0;
        end
        rfoptions{3} = 'microscope';
        handles.refFrameVal = rfoptions{idx};
    end
  
end
  
function plotDelCB(hObj,event)
    
  hs = handles;
  
  % check whether or not there is any data
  if length(hs.data)==1 && isempty(hs.data{1})
      errorbox('Please add some experiments before attempting to plot.');
      return
  end
  
  % get the legend
  legend = getLegend;
  
  % get which measure has been requested
  meas = hObj.String(18);
  switch meas
      case 'x'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','deltaX',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case 'y'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','deltaY',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case 'z'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','deltaZ',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '1'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta1D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '2'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta2D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '3'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta3D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
  end
  set(gca,'FontSize',12);
  
  % enable save and export buttons
  handles.savebtn.Enable = 'on';
  handles.exportbtn.Enable = 'on';
    
end

function plotKKCB(hObj,event)
    
  hs = handles;
  
  % check whether or not there is any data
  if length(hs.data)==1 && isempty(hs.data{1})
      errorbox('Please add some experiments before attempting to plot.');
  end
  
  % get the legend  
  legend = getLegend;
  
  % get which measure has been requested
  meas = hObj.String(13);
  switch meas
      case 'x'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','sisSepX',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case 'y'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','sisSepY',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case 'z'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','sisSepZ',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '1'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','sisSep1D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '2'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','sisSep2D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '3'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','sisSep3D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
  end
  set(gca,'FontSize',12);
  
  % enable save and export buttons
  handles.savebtn.Enable = 'on';
  handles.exportbtn.Enable = 'on';
    
end

function plotSwivCB(hObj,event)
    
  hs = handles;
  
  % check whether or not there is any data
  if length(hs.data)==1 && isempty(hs.data{1})
      errorbox('Please add some experiments before attempting to plot.');
  end
  
  % get the legend 
  legend = getLegend;
  
  % get which measure has been requested
  meas = hObj.String(18);
  switch meas
      case 'y'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','swivelY',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case 'z'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','swivelZ',...
              'depthFilter',hs.filtering.Value,'legend',legend);
      case '3'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','swivel3D',...
              'depthFilter',hs.filtering.Value,'legend',legend);
  end
  set(gca,'FontSize',12);
  
  % enable save and export buttons
  handles.savebtn.Enable = 'on';
  handles.exportbtn.Enable = 'on';
    
end

function plot2DcorrCB(hObj,event)

  hs = handles;
  
  % check whether or not there is any data
  if isempty(hs.data{1})
      errorbox('Please add some experiments before attempting to plot.');
  end
  
  % get which measure has been requested
  meas = hObj.String(18);
  switch meas
      case '<'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta2DsisSep2D',...
              'depthFilter',hs.filtering.Value);
      case 't'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta2DswivelY',...
              'depthFilter',hs.filtering.Value);
  end
  set(gca,'FontSize',12);
  
  % enable save and export buttons
  handles.savebtn.Enable = 'on';
  handles.exportbtn.Enable = 'on';

end

function plot3DcorrCB(hObj,event)

  hs = handles;
  
  % check whether or not there is any data
  if isempty(hs.data{1})
      errorbox('Please add some experiments before attempting to plot.');
  end
  
  % get which measure has been requested
  meas = hObj.String(18);
  switch meas
      case '<'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta3DsisSep3D',...
              'depthFilter',hs.filtering.Value);
      case 't'
          basicPlots_old(hs.data,'refFrame',hs.refFrameVal,'stat','delta3Dswivel3D',...
              'depthFilter',hs.filtering.Value);
  end
  set(gca,'FontSize',12);
  
  % enable save and export buttons
  handles.savebtn.Enable = 'on';
  handles.exportbtn.Enable = 'on';

end

function saveCB(hObj,event)
  
  % ask the user to select a directory
  [filename,savepath] = uiputfile('*.eps','Save figure','figure.eps');
  
  % open a new figure, and copy all figPanel objects
  h2 = figure;
  copyobj(get(handles.figPanel,'children'),h2);
  print(h2,fullfile(savepath,filename),'-depsc');
  close(h2);
  
end

function exportCB(hObj,event)
  
    % open a new figure, and copy all figPanel objects
    h2 = figure;
    copyobj(get(handles.figPanel,'children'),h2);
  
end

%% Other functions.

function legend = getLegend
  
  % loop over experiments to get legends
  for iExpt = 1:handles.nExpts
    legend{iExpt} = handles.exptLabel{iExpt}.String;
    % if no label given, replace with expt #
    if isempty(legend{iExpt})
        legend{iExpt} = ['Expt ' num2str(iExpt)];
    end
  end
  
end

function data = loadData(filepath)
  % get file
  data = load(filepath);
  nmes = fieldnames(data);
  if length(nmes)>1
    errorbox('Cannot load data: file contains more than one variable.')
    data = [];
  else
    data = getfield(data,nmes{1});
  end
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

end % kitIntensityAnalysis
