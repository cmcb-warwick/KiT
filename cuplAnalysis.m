function cuplAnalysis(jobset,ch)
% CUPLANALYSIS Produces a GUI using which the user can produce multiple
% plots.
%
%   cuplAnalysis(jobset,channel)
%
%  jobset: Struct containing tracking job setup options, or cell array of
%  multiple jobsets.
%
%  channel: Channel to plot. Defaults to coordinate system channel.
%
% Created by: Jonathan W. Armond
% Edited by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

if ~iscell(jobset)
  jobset = {jobset};
end
njs = length(jobset);

if nargin<2
  ch = jobset{1}.options.coordSystemChannel;
end

% Upgrade jobset if required.
for i=1:njs
  if ~isfield(jobset{i}.options,'jobsetVersion') || ...
    jobset{i}.options.jobsetVersion < kitVersion(2)
    jobset{i} = kitJobset(jobset{i});
  end
end
% Number of ROIs in each job.
for i=1:njs
  nroi(i) = length(jobset{i}.ROI);
end

% Setup GUI.
roiChanged = 1;
analysisLoaded = 0;
handles = createControls(jobset);
ROICB();
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls(jobset)
  
  % Define font sizes (default is 12)
  largefont = 16;
  medfont = 14;
  tinyfont = 10;
    
  % Define some standard measures.
  dx = 2;
  ddx = 1;
  h = 2;
  lh = 1.5; % label height
    
  % Create figure.
  colw = 40;
  figw = 2*colw+3*dx;
  figh = 42;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = 'CupL analysis';
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';
  movegui(hs.fig,'center');
  
  % Define component positions.
  lw = 0.85*colw; % label width
  ex = lw + ddx; % edit box start
  ew = colw - lw; % edit box width
  mw = 0.25*colw; % drop menu width
  
  % Define x and y positions.
  x = dx;
  y = figh-(ddx+lh);
  
  t=label(hs.fig,'Jobset:',[x y colw lh],largefont);
  t.FontWeight='bold';
  y=y-lh;
  [~,name] = fileparts(jobset{1}.filename);
  name = strshorten(name,60);
  label(hs.fig,name,[x y figw lh],medfont);
  y=y-h;
  label(hs.fig,['Reading channel ' num2str(ch) ' data'],[x y colw lh]);
  y=y-h;
  
  savedy = y; % will come back here later
  
  % Per-movie column.
  t=label(hs.fig,'Per-movie analysis',[x y colw lh],medfont);
  t.FontWeight='bold'; t.HorizontalAlignment = 'Center';

  % Get ROI names.
  maxMovLen = 32;
  ROIString = {};
  idx = 1;
  for j=1:njs
    for i=1:length(jobset{j}.ROI)
      ROIString{idx} = sprintf('%d: %s',idx,strshorten(jobset{j}.ROI(i).movie,maxMovLen));
      idx = idx+1;
    end
  end
  y=y-h;
  label(hs.fig,'Selected movie:',[x y colw lh]);
  y=y-h;
  hs.ROI = popup(hs.fig,ROIString,[x y colw h],@ROICB);

  y=y-h;
  label(hs.fig,'Min length %',[x y colw lh]);
  hs.minLength = editbox(hs.fig,'75',[ex y ew lh]);
  y=y-h;
  hs.sisters = button(hs.fig,'Sister trajectories',[x y colw h],@sistersCB);
  y=y-h;
  hs.tracks = button(hs.fig,'Track overlay',[x y colw h],@tracksCB);
  y=y-h;
  hs.plate = button(hs.fig,'Metaphase plate',[x y colw h],@plateCB);
  y=y-h;
  t=label(hs.fig,'Channel',[x y mw lh]); t.HorizontalAlignment = 'Center';
  y=y-h;
  hs.intensityCh = popup(hs.fig,'1',[x y mw h]);
  hs.intensity = button(hs.fig,'Intensity/position trajectories',[x+(mw+ddx) y colw-(mw+ddx) h],@intensityCB);
  y=y-h;
  t=label(hs.fig,'Timepoint',[x y mw lh]); t.HorizontalAlignment = 'Center';
  t=label(hs.fig,'Sister pair',[x+(mw+ddx) y mw lh]); t.HorizontalAlignment = 'Center';
  y=y-h;
  % SHOW IMAGE
  hs.imageTp = popup(hs.fig,'1',[x y mw h]);
  t=label(hs.fig,'-',[x+(mw+ddx) y mw lh]); t.HorizontalAlignment = 'Center';
  hs.image = button(hs.fig,'Show full image',[x+2*(mw+ddx) y colw-2*(mw+ddx) h],@imageCB);
  y=y-h;
  % SHOW ALL SPOTS
  hs.allSpotsTp = popup(hs.fig,'1',[x y mw h]);
  t=label(hs.fig,'-',[x+(mw+ddx) y mw lh]); t.HorizontalAlignment = 'Center';
  hs.allSpots = button(hs.fig,'Show all spots',[x+2*(mw+ddx) y colw-2*(mw+ddx) h],@allSpotsCB);
  y=y-h;
  % SHOW ALL SISTERS
  hs.allSisTp = popup(hs.fig,'1',[x y mw h]);
  t=label(hs.fig,'-',[x+(mw+ddx) y mw lh]); t.HorizontalAlignment = 'Center';
  hs.allSis = button(hs.fig,'Show all sisters',[x+2*(mw+ddx) y colw-2*(mw+ddx) h],@allSisCB);
  y=y-h;
  % SHOW SISTER PAIR
  hs.sisPairTp = popup(hs.fig,'1',[x y mw h]);
  hs.sisPairSis = popup(hs.fig,'1',[x+(mw+ddx) y mw h]);
  hs.sisPair = button(hs.fig,'Show sister pair',[x+2*(mw+ddx) y colw-2*(mw+ddx) h],@sisPairCB);
  y=y-h;
  % DRAGON TAILS
  hs.dragTailTp = popup(hs.fig,'1',[x y mw h]);
  hs.dragTailSis = popup(hs.fig,'1',[x+(mw+ddx) y mw h]);
  hs.dragTail = button(hs.fig,'Show dragon tails',[x+2*(mw+ddx) y colw-2*(mw+ddx) h],@dragTailCB);
  y=y-h;
  % SHOW IMAGE
  t=label(hs.fig,'-',[x y mw lh]); t.HorizontalAlignment = 'Center';
  hs.kymoSis = popup(hs.fig,'1',[x+(mw+ddx) y mw h]);
  hs.kymo = button(hs.fig,'Show kymograph',[x+2*(mw+ddx) y colw-2*(mw+ddx) h],@kymoCB);
  
  % Re-define all positions for second column.
  x = colw+2*dx;
  y=savedy;
  ex = x + ex; % edit box start
  
  % Experiment column.
  t=label(hs.fig,'Experiment analysis',[x y colw lh],medfont);
  t.FontWeight='bold'; t.HorizontalAlignment = 'Center';
  y=y-h;
  hs.diags = button(hs.fig,'Print diagnostics',[x y colw h],@diagsCB);
  y=y-h;
  hs.autocorr = button(hs.fig,'Autocorrelation',[x y colw h],@autocorrCB);
  y=y-h;
  hs.dists = button(hs.fig,'Intersister distance',[x y colw h],@distsCB);
  y=y-h;
  hs.speeds = button(hs.fig,'Speeds',[x y colw h],@speedsCB);
  y=y-h;
  hs.msd = button(hs.fig,'Mean-squared displacement',[x y colw h],@msdCB);
  y=y-h;
  hs.spindle = button(hs.fig,'Spindle length',[x y colw h],@spindleCB);
  y=y-2*h;
  hs.save = button(hs.fig,'Save analysis',[x y colw h],@saveCB);
  y=y-h;
  hs.status = label(hs.fig,'',[x y colw lh]);
  hs.status.ForegroundColor = [1 0 0];
  
  % Close button.
  btnw = 10;
  y=1; x = figw-(btnw+dx);
  hs.close = button(hs.fig,'Close',[x y btnw h],@closeCB);
  
end

function ROICB(hObj,event)
  roiChanged = 1;
  job = loadActiveJob();
  % WANT TO CONSTRUCT TIMEPOINT, SISTER AND CHANNEL NUMBERS
  chans = {1:4};
  handles.intensityCh.String = chans;
  time = {1:job.metadata.nFrames};
  handles.imageTp.String = time;
  handles.allSpotsTp.String = time;
  handles.allSisTp.String = time;
  handles.sisPairTp.String = time;
  handles.dragTailTp.String = time;
  sis = {1:length(job.dataStruct{ch}.sisterList)};
  handles.sisPairSis.String = sis;
  handles.dragTailSis.String = ['All';sis];
  handles.kymoSis.String = sis;
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
  signalCh = str2double(handles.intensityCh.String(handles.intensityCh.Value));
  if ~job.options.intensity.execute(signalCh)
      errormsg = ['Channel ' num2str(signalCh) ' contains no intensity information'];
      errorbox(errormsg);
      return
  end
  minLen = str2double(handles.minLength.String)/100;
  kitPlotSisterSignal(job,signalCh,'channel',ch,'minLength',minLen);
end

function imageCB(hObj,event)
  job = loadActiveJob();
  imageTp = str2double(handles.imageTp.String(handles.imageTp.Value));
  chans = 1:job.metadata.nChannels;
  
  chanLabs = 'bgrm';
  wvs = job.metadata.wavelength*1000;
  for iwv = 1:length(wvs)
    if isnan(wvs(iwv))
      chanOrder(iwv)='n';
    elseif wvs(iwv)<490
      chanOrder(iwv)='b';
    elseif wvs(iwv)<580
      chanOrder(iwv)='g';
    elseif wvs(iwv)<670
      chanOrder(iwv)='r';
    else
      chanOrder(iwv)='m';
    end
  end
  
  updateStatus('Generating image');
  kitShowImage(job,'imageChans',chans,'timePoint',imageTp,'chanOrder',chanOrder);
  updateStatus('');
end

function allSpotsCB(hObj,event)
  job = loadActiveJob();
  allSpotsTp = str2double(handles.allSpotsTp.String(handles.allSpotsTp.Value));
  chan = job.options.coordSystemChannel;
  
  updateStatus('Generating all spots');
  kitShowSpots(job,'channel',chan,'timePoint',allSpotsTp);
  updateStatus('');
end

function allSisCB(hObj,event)
  job = loadActiveJob();
  allSisTp = str2double(handles.allSisTp.String(handles.allSisTp.Value));
  chan = job.options.coordSystemChannel;
  
  updateStatus('Generating all sisters');
  kitShowSisters(job,'channel',chan,'timePoint',allSisTp);
  updateStatus('');
end

function sisPairCB(hObj,event)
  job = loadActiveJob();
  sisPairTp = str2double(handles.sisPairTp.String(handles.sisPairTp.Value));
  sisPairSis = str2double(handles.sisPairSis.String(handles.sisPairSis.Value));
  chan = job.options.coordSystemChannel;
  
  updateStatus('Generating image of sister pair');
  kitShowSisterPair(job,'channel',chan,'timePoint',sisPairTp,'sisterPair',sisPairSis);
  updateStatus('');
end

function dragTailCB(hObj,event)
  job = loadActiveJob();
  dragTailTp = str2double(handles.dragTailTp.String(handles.dragTailTp.Value));
  dragTailSis = str2double(handles.dragTailSis.String(handles.dragTailSis.Value));
  chan = job.options.coordSystemChannel;
  
  if isnan(dragTailSis)
    trks = [];
  else
    trks = job.dataStruct{chan}.sisterList(1).trackPairs(dragTailSis,1:2);
  end
  
  updateStatus('Generating image of dragon tails');
  kitShowDragonTails(job,'channel',chan,'timePoint',dragTailTp,'tracks',trks);
  updateStatus('');
end

function kymoCB(hObj,event)
  job = loadActiveJob();
  kymoSis = str2double(handles.kymoSis.String(handles.kymoSis.Value));
  
  wvs = job.metadata.wavelength*1000;
  for iwv = 4:-1:1
    if isnan(wvs(iwv))
      continue
    elseif wvs(iwv)<490
      chanMap(iwv)=3;
    elseif wvs(iwv)<580
      chanMap(iwv)=2;
    elseif wvs(iwv)<670
      chanMap(iwv)=1;
    end
  end
  
  updateStatus('Generating kymograph');
  kitMakeKymograph(job,kymoSis,'channelMap',chanMap);
  updateStatus('');
end

function diagsCB(hObj,event)
  for j=1:njs
    s = ['Generating jobset diagnostics: jobset ' num2str(j)];
    updateStatus(s);
    kitJobsetDiagnostics(jobset{j},ch);
  end
  updateStatus('');
end

function autocorrCB(hObj,event)
  analysis = loadAnalysis();
  updateStatus('Computing autocorrelations');
  analysis = cuplAutocorrelation(analysis);
  updateStatus('');
  autocorrs = analysis.autocorrs;
  [~,locs] = findpeaks(-autocorrs.sisters.m_dx,autocorrs.t);  
  fprintf('Minima in the autocorrelation plot are at %f suggesting the period is %f s\n',locs(1),2*locs(1));
  cuplPlotCorrelation(autocorrs.t,autocorrs.sisters.m_dx,autocorrs.sisters.e_dx,...
    'PlotTitle','Mean autocorrelation of sister pair centres, dx');
  cuplPlotCorrelation(autocorrs.t,autocorrs.sisters.cm_dx,autocorrs.sisters.ce_dx, ...
    'PlotTitle','Mean autocorrelation by cell of sister pair centres, dx');
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

function closeCB(hObj,event)
  uiresume(gcf);
end

function job = loadActiveJob()
  persistent jobP
  if roiChanged
    idx = handles.ROI.Value;
    % Figure out which jobset is selected.
    sroi = cumsum(nroi);
    jsidx = find(sroi>=idx,1);
    if jsidx>1
      idx = idx-sroi(jsidx-1);
    end
    updateStatus('Loading ROI');
    jobP = kitLoadJob(jobset{jsidx},idx);
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

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end
