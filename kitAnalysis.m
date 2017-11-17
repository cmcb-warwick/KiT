function kitAnalysis(jobset,ch)
% KITBASICPLOTS Produces a set of basic plots for a job
%
%   kitAnalysis(jobset,channel)
%
%  jobset: Struct containing tracking job setup options.
%
%  channel: Channel to plot. Defaults to 1.
%
% Copyright (c) 2015 Jonathan W. Armond

if nargin<2
  ch = 1;
end

% Upgrade jobset if required.
if ~isfield(jobset.options,'jobsetVersion') || ...
    jobset.options.jobsetVersion < kitVersion(2)
  jobset = kitJobset(jobset);
end

% Setup GUI.
roiChanged = 1;
analysisLoaded = 0;
handles = createControls(jobset);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls(jobset)
  w = 30;
  h = 54;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 w+2 h]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiT ' kitVersion(1)];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';
  movegui(hs.fig,'center');
  figpos = hs.fig.Position;

  x = 1;
  logoh = 8;
  y = figpos(4)-logoh-1;
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[x y w logoh]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kitlogo.png'),pos([4 3])));

  h = 2;
  y = y-h;
  lh = 1.5; % label height
  lw = 0.75*w; % label width
  ex = lw + 1; % edit box start
  ew = 0.25*w - 1; % edit box width

  [~,name] = fileparts(jobset.filename);
  t=label(hs.fig,['Jobset: ' name],[x y w lh]);
  t.FontWeight='bold';
  y=y-h;
  label(hs.fig,['Reading channel ' num2str(ch) ' data'],[x y w lh],10);
  y=y-h;
  t=label(hs.fig,'Per ROI analysis',[x y w lh]);
  t.FontWeight='bold';

  % Get ROI names.
  maxMovLen = 32;
  ROIString = {};
  for i=1:length(jobset.ROI)
    ROIString{i} = sprintf('%d: %s [%s]',i,strshorten(jobset.ROI(i).movie,maxMovLen),num2str(round(jobset.ROI(i).crop),'%d '));
  end
  y=y-h;
  label(hs.fig,'ROIs:',[x y w lh]);
  y=y-h;
  hs.ROI = popup(hs.fig,ROIString,[x y w h],@ROICB);

  y=y-h;
  label(hs.fig,'Min length %',[x y w lh]);
  hs.minLength = editbox(hs.fig,'75',[ex y ew h]);
  y=y-h;
  hs.sisters = button(hs.fig,'Sister trajectories',[x y w h],@sistersCB);
  y=y-h;
  hs.tracks = button(hs.fig,'Track overlay',[x y w h],@tracksCB);
  y=y-h;
  hs.plate = button(hs.fig,'Metaphase plate',[x y w h],@plateCB);
  y=y-h;
  label(hs.fig,'Intensity channel',[x y w lh]);
  hs.intensityCh = editbox(hs.fig,'1',[ex y ew h]);
  y=y-h;
  hs.intensity = button(hs.fig,'Intensity/position trajectories',[x y w h],@intensityCB);

  y=y-2*h;
  t=label(hs.fig,'Experiment analysis',[x y w lh]);
  t.FontWeight='bold';
  y=y-h;
  hs.diags = button(hs.fig,'Print diagnostics',[x y w h],@diagsCB);
  y=y-h;
  hs.autocorr = button(hs.fig,'Autocorrelation',[x y w h],@autocorrCB);
  y=y-h;
  hs.dists = button(hs.fig,'Intersister distance',[x y w h],@distsCB);
  y=y-h;
  hs.speeds = button(hs.fig,'Speeds',[x y w h],@speedsCB);
  y=y-h;
  hs.msd = button(hs.fig,'Mean-squared displacement',[x y w h],@msdCB);
  y=y-h;
  hs.spindle = button(hs.fig,'Spindle length',[x y w h],@spindleCB);
  y=y-2*h;
  hs.save = button(hs.fig,'Save analysis',[x y w h],@saveCB);
  y=y-h;
  hs.status = label(hs.fig,'',[x y w lh]);
end

function ROICB(hObj,event)
  roiChanged = 1;
end

function sistersCB(hObj,event)
  job = loadActiveJob();
  minLen = str2double(handles.minLength.String)/100;
  kitPlotSisters(job,'minLength',minLen,'channel',ch);
end

function tracksCB(hObj,event)
  job = loadActiveJob();
  minLen = str2double(handles.minLength.String)/100;
  kitPlotTracks(job,'minLength',minLen,'channel',ch);
end

function plateCB(hObj,event)
  job = loadActiveJob();
  kitPlateCheck(job,'channel',ch);
end

function intensityCB(hObj,event)
  job = loadActiveJob();
  signalCh = str2double(handles.intensityCh.String);
  minLen = str2double(handles.minLength.String)/100;
  kitPlotSisterSignal(job,signalCh,'channel',ch,'minLength',minLen);
end

function diagsCB(hObj,event)
  kitJobsetDiagnostics(jobset,ch);
end

function autocorrCB(hObj,event)
  analysis = loadAnalysis();
  updateStatus('Computing autocorrelations');
  analysis = cuplAutocorrelation(analysis);
  updateStatus('');
  autocorrs = analysis.autocorrs;
  cuplPlotCorrelation(autocorrs.t,autocorrs.sisters.m_dx,autocorrs.sisters.e_dx,...
    'PlotTitle','Mean autocorrelation of sister pair centres, dx');
end

function distsCB(hObj,event)
  analysis = loadAnalysis();
  updateStatus('Computing sister distances');
  analysis = cuplDistances(analysis);
  updateStatus('');
  cuplPlotHistogram(analysis.distances.sisters.sister.d,'PlotTitle','Distance between sisters',...
                    'PlotXLabel','Distance (um)','plotYLabel','Density','NumBins',20);
end

function speedsCB(hObj,event)
  analysis = loadAnalysis();
  updateStatus('Computing speeds');
  analysis = cuplSpeeds(analysis);
  updateStatus('');
  speeds = analysis.speeds.sisters.dx*60; % Convert to um/min.
  cuplPlotHistogram(speeds,'plotYLabel','Density',...
                    'plotXLabel','Normal speed (um/min)','PlotTitle','Mean speeds')
end

function msdCB(hObj,event)
  analysis = loadAnalysis();
  updateStatus('Computing mean square displacement');
  analysis = cuplMsd(analysis);
  updateStatus('');
  msd = analysis.msd.time;
  cuplPlotBins(analysis.time,msd.m_d,msd.e_d,'PlotTitle','Mean squared displacement',...
               'PlotXLabel','Time (s)','PlotYLabel','Mean square distance (\mum^2)');
end

function spindleCB(hObj,event)
  analysis = loadAnalysis();
  updateStatus('Computing spindle lengths');
  analysis = cuplSpindleLength(analysis);
  updateStatus('');
  len = analysis.spindle.length;
  cuplPlotHistogram(len,'PlotTitle','Pole-pole distance',...
                    'PlotXLabel','Distance (um)','plotYLabel','Density','NumBins',20);
end

function saveCB(hObj,event)
  [filename,pathname] = uiputfile('*.mat','Save file name for analysis','analysis.mat');
  if ~(isequal(filename,0) || isequal(pathname,0))
    % Run analyses.
    analysis = loadAnalysis();
    updateStatus('Computing autocorrelations');
    analysis = cuplAutocorrelation(analysis);
    updateStatus('Computing speeds');
    analysis = cuplSpeeds(analysis);
    updateStatus('Computing mean square displacement');
    analysis = cuplMsd(analysis);
    updateStatus('Computing sister distances');
    analysis = cuplDistances(analysis);
    updateStatus('Computing spindle length');
    analysis = cuplSpindleLength(analysis);
    analysis = cuplStripData(analysis);

    updateStatus('Saving');
    save(fullfile(pathname,filename),'analysis');
    updateStatus('Saved');
  end
end

function job = loadActiveJob()
  persistent jobP
  if roiChanged
    idx = handles.ROI.Value;
    updateStatus('Loading ROI');
    jobP = kitLoadJob(jobset,idx);
    updateStatus('');
    roiChanged = 0;
  end
  job = jobP;
end

function analysis = loadAnalysis()
  persistent analysisP;
  if ~analysisLoaded || isempty(analysisP)
    analysisP.options.percentNan = 1 - str2double(handles.minLength.String)/100;
    analysisP.options.minSistersPerCell = 1;
    analysisP.options.byDistNumBins = 12;
    analysisP.options.byDistMaxWidth = 12;
    analysisP.options.neighbourThreshold = 3;
    analysisP.options.monotelicAngle = 8;
    analysisP.options.trackChannel = 1;
    analysisP.options.channel = 1;
    analysisP.options.correctDrift = 0;
    analysisP.options.keepAllData = 1;
    analysisP.options.doTracks = 1;
    analysisP.options.poleCutoff = 4;
    analysisP.stages = {};

    updateStatus('Loading data');
    analysisP = cuplLoadData(analysisP,jobset,ch);
    updateStatus('Preprocessing data');
    analysisP = cuplPreprocess(analysisP);
    updateStatus('');
    analysisLoaded = 1;
  end
  analysis = analysisP;
end

function updateStatus(s)
  handles.status.String = s;
  drawnow;
end

end % kitAnalysis

%% LOCAL FUNCTIONS

function m=expand(m,rows)
% Expand m to have at least rows.
  sizeDiff = nFrames - size(m,1);
if ~isempty(m) && sizeDiff > 0
  m = [m; nan(sizeDiff,size(m,2))];
end
end
