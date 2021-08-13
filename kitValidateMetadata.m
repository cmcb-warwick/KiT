function [jobset,applyAll] = kitValidateMetadata(jobset,iMov)
% KITVALIDATEMETADATA Allows user to check and correct automatically-obtained
% metadata.
%
% Copyright (c) 2017 C. A. Smith

% get movie filepath
movieDirectory = jobset.movieDirectory;
movie = jobset.ROI(iMov).movie;
movieFilePath = [movieDirectory '/' movie];

%% GET METADATA FROM MOVIE

if isfield(jobset,'metadata') && length(jobset.metadata)>=iMov && jobset.metadata{iMov}.validated
    
    % if there is already validated metadata, get it
    md = jobset.metadata{iMov};
    
else
    % otherwise, get metadata from the movie
    kitLog('Opening movie: %s', movieFilePath);
    if ~exist(movieFilePath,'file')
      error('Could not find file: %s',movieFilePath);
    end

    addpath bfmatlab;
    bfCheckJavaPath(1);

    reader = bfGetReader(movieFilePath);

    % Read basic image metadata.
    md.nFrames = reader.getSizeT();
    md.nChannels = reader.getSizeC();
    md.nPlanes = reader.getImageCount();
    pixelType = reader.getPixelType();
    md.nBytesPerPixel = loci.formats.FormatTools.getBytesPerPixel(pixelType);
    md.isSigned = loci.formats.FormatTools.isSigned(pixelType);
    md.isFloatingPoint = loci.formats.FormatTools.isFloatingPoint(pixelType);
    md.isLittleEndian = reader.isLittleEndian();
    switch md.nBytesPerPixel
      case 1
        if md.isSigned, md.dataType = 'int8'; else md.dataType = 'uint8'; end
      case 2
        if md.isSigned, md.dataType = 'int16'; else md.dataType = 'uint16'; end
      case 4
        md.dataType = 'single';
      case 8
        md.dataType = 'double';
      otherwise
        error('Only 8- or 16-bit integer, or 32- or 64-bit floating point data supported');
    end
    md.frameSize = [reader.getSizeX() reader.getSizeY() reader.getSizeZ()];

    % Read additional metadata
    metaTable = reader.getMetadataStore();
    try
      numWvs = metaTable.getChannelCount(0);
    catch
      numWvs = 3;
    end
    md.wavelength = [525 615 705]/1000; % Default assumes eGFP, mCherry, far-red. FIXME Ask user.
    for iChan=1:numWvs
      try
        chWv = metaTable.getChannelEmissionWavelength(0,iChan-1);
        try
          md.wavelength(iChan) = chWv.value(ome.units.UNITS.MICROM).doubleValue();
        catch
          md.wavelength(iChan) = chWv.getValue();
        end
      catch
        continue
      end
    end
    % Timepoints per plane
    nZPlanes = md.frameSize(3);
    nTimepoints = md.nFrames;
    counter = 0;
    defDt = 2; % Default assume every 2 sec. FIXME Ask user.
    for iTime=1:nTimepoints
      defT = (iTime-1)*defDt;
      for iPlane=1:nZPlanes
        try
          md.frameTime(iPlane,iTime) = metaTable.getPlaneDeltaT(0, counter).doubleValue();
        catch
	  try
	    md.frameTime(iPlane,iTime) = metaTable.getPlaneDeltaT(0, counter).value.double/1000; % in this format provided in milliseconds
	  catch
            % Use default, if missing metadata.
            md.frameTime(iPlane,iTime) = defT;
	  end
        end
        counter = counter+1;
      end
    end
    dT = diff(md.frameTime(1,:));
    % Is 3D image?
    md.is3D = nZPlanes > 1;
    % Numerical aperture
    try
      md.na = metaTable.getObjectiveLensNA(0,0).doubleValue;
    catch
      md.na = 1.4; % Assume default.
    end
    % Physical pixel size
    try
      try
        if md.is3D
          md.pixelSize = [
            metaTable.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM).doubleValue(),...
            metaTable.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM).doubleValue(),...
            metaTable.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM).doubleValue()];
        else
          md.pixelSize = [
            metaTable.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM).doubleValue(),...
            metaTable.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM).doubleValue(),...
            1];
        end
      catch
        if md.is3D
          md.pixelSize = [
            metaTable.getPixelsPhysicalSizeX(0).getValue(),...
            metaTable.getPixelsPhysicalSizeY(0).getValue(),...
            metaTable.getPixelsPhysicalSizeZ(0).getValue()];
        else
          md.pixelSize = [
            metaTable.getPixelsPhysicalSizeX(0).getValue(),...
            metaTable.getPixelsPhysicalSizeY(0).getValue(),...
            1];
        end
      end
    catch
      if md.is3D 
        md.pixelSize = [0.0645 0.0645 0.2];
      else
        md.pixelSize = [0.0645 0.0645 1];
      end
    end
    md.validated = 0;
end

% back up metadata in case
jobset.metadata{iMov} = md;

%% GUI

% Setup GUI.
handles = createControls();
updateControls(md);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls()
  
  figw = 50;
  figh = 29;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = 'Check metadata';
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  medfont = 14;
  smallfont = 12;
  tinyfont = 10;
  
  % Set some standard positions and distances
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5;
  ddx = 1;
  toplabely = figh-lh; %top-most point
  
  % Editable jobset name at top of the screen
  colw = figw-2*dx;
  hs.movieNameText = label(hs.fig,'Movie file name',[dx toplabely colw h],medfont);
  hs.movieNameText.FontWeight = 'bold';
  y = toplabely-lh;
  hs.movieName = label(hs.fig,[],[dx+ddx y colw lh],smallfont);
  
  %% List of metadata for validation
  
  % dimensions panel
  panx = 2*dx; panw = colw-panx;
  panh = 6.25*h; pany = y-(panh+0.5*h);
  hs.dimensionsPanel = uipanel(hs.fig,'Units','characters','Position',[panx pany panw panh],'FontSize',tinyfont,'Title','Movie dimensions');
  p = hs.dimensionsPanel;
  
  % Frame size, z-slices and t-points
  labw = 25;
  editx = dx+(labw+ddx);
  editw = panw-(editx+dx);
  edity = panh-(2*h);
  hs.nXYpixelsText = label(p,'Frame size: [x,y]',[dx edity labw h],tinyfont);
  hs.nXYpixels = label(p,[],[editx edity editw h],tinyfont);
  edity = edity-h;
  hs.is3D = checkbox(p,'3D movies?',[dx edity labw h],@is3DCB,tinyfont);
  edity = edity-h;
  hs.nPlanesText = label(p,'z-slices',[dx edity labw h],tinyfont);
  hs.nPlanes = editbox(p,[],[editx edity editw h],tinyfont);
  edity = edity-h;
  hs.nFramesText = label(p,'Time points',[dx edity labw h],tinyfont);
  hs.nFrames = editbox(p,[],[editx edity editw h],tinyfont);
  edity = edity-h;
  
  % Channels
  nChans = 3;
  labw = 10;
  editx = dx+(labw+ddx);
  editw = (panw-(editx+dx+(nChans-1)*ddx))/nChans;
  hs.nChannelsText = label(p,'Channels',[dx edity labw h],tinyfont);
  for iChan=1:nChans
    hs.channelNum{iChan} = uicontrol('Parent',p,'Units','characters','Style','radio','String',num2str(iChan),'Position',[editx edity editw h],'Callback',@nChannelsCB);
    editx = editx+(editw+ddx);
  end
  hs.channelNum{1}.Value = 1;
  hs.nChannels = 1;
  
  % Wavelengths
  y = y-(panh+3*h);
  labw = 15;
  editx = dx+(labw+ddx);
  editw = (colw-(editx+(nChans-1)*ddx))/nChans;
  hs.wavelengthText = label(hs.fig,'Wavelength (nm)',[dx y labw h],tinyfont);
  for iChan=1:nChans
    hs.waveChanText{iChan} = label(hs.fig,['Ch.' num2str(iChan)],[editx y+0.75*h editw h],tinyfont);
    hs.waveChanText{iChan}.HorizontalAlignment = 'center';
    hs.waveChanText{iChan}.FontWeight = 'bold';
    hs.wavelength{iChan} = editbox(hs.fig,[],[editx y editw h],tinyfont);
    editx = editx+(editw+ddx);
  end
  
  % Pixel sizes
  nCoords = 3;
  y = y-2*h;
  labw = 15;
  editx = dx+(labw+ddx);
  editw = (colw-(editx+(nCoords-1)*ddx))/nCoords;
  coordLabel = ['x','y','z'];
  hs.pixelSizeText = label(hs.fig,'Pixel size (nm)',[dx y labw h],tinyfont);
  for iCoord=1:nCoords
    hs.pixCoordText{iCoord} = label(hs.fig,coordLabel(iCoord),[editx y+0.75*h editw h],tinyfont);
    hs.pixCoordText{iCoord}.HorizontalAlignment = 'center';
    hs.pixCoordText{iCoord}.FontWeight = 'bold';
    hs.pixelSize{iCoord} = editbox(hs.fig,[],[editx y editw h],tinyfont);
    editx = editx+(editw+ddx);
  end
  
  % Time lapse and NA
  y = y-lh;
  labw = 30;
  editx = dx+(labw+ddx);
  editw = colw-(editx);
  hs.timeLapseText = label(hs.fig,'Time lapse, dt (s)',[dx y labw h],tinyfont);
  hs.timeLapse = editbox(hs.fig,[],[editx y editw h],tinyfont);
  y = y-h;
  hs.numAperText = label(hs.fig,'Numerical aperture',[dx y labw h],tinyfont);
  hs.numAper = editbox(hs.fig,[],[editx y editw h],tinyfont);
  
  % Apply to all
  y = y-h;
  hs.applyAll = checkbox(hs.fig,'Apply to all movies?',[dx y labw h],[],tinyfont);
  
  % Validate button
  btnw = 10; btnh = 2;
  y = 1.25;
  btnx = figw-(btnw+dx);
  hs.validate = button(hs.fig,'Validate',[btnx y btnw btnh],@saveCB);
  
  movegui(hs.fig,'center');
  
end

%% Update control status based on contents of metadata.
function updateControls(md)
  hs = handles;
  
  idx = find(jobset.ROI(1).movie=='/',1,'last');
  if isempty(idx)
    idx=0;
  end
  hs.movieName.String = jobset.ROI(iMov).movie(idx+1:end);
  hs.nXYpixels.String = [num2str(md.frameSize(1)) ' x ' num2str(md.frameSize(2))];
  hs.is3D.Value = md.is3D;
  hs.nPlanes.String = num2str(md.nPlanes/(md.nChannels*md.nFrames));
  hs.nFrames.String = num2str(md.nFrames);
  hs.nChannels = md.nChannels;
  for iChan=1:length(md.wavelength)
    hs.channelNum{iChan}.Value = (iChan<=hs.nChannels);
    hs.wavelength{iChan}.String = num2str(md.wavelength(iChan)*1000);
  end
  for iCoord=1:3
    hs.pixelSize{iCoord}.String = num2str(md.pixelSize(iCoord)*1000);
  end
  if str2double(hs.nFrames.String) == 1
    hs.timeLapse.String = 'N/A';
  else
    hs.timeLapse.String = num2str(md.frameTime(1,2));
  end
  hs.numAper.String = num2str(md.na);
  
  handles = hs;
  
  nChannelsCB();
  
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

function nChannelsCB(hObj,event)
  if exist('hObj','var')
    nChans = str2double(hObj.String);
    handles.nChannels = nChans;
    for notChan = setdiff(1:3,nChans)
      handles.channelNum{notChan}.Value = 0;
    end
  end
  for iChan = 1:3
    if iChan > handles.nChannels
      handles.wavelength{iChan}.Enable = 'off';
      handles.waveChanText{iChan}.Enable = 'off';
    else
      handles.wavelength{iChan}.Enable = 'on';
      handles.waveChanText{iChan}.Enable = 'on';
    end
  end
end

function is3DCB(hObj,event)
  if handles.is3D.Value
    handles.nPlanes.Enable = 'on';
    handles.pixCoordText{3}.Enable = 'on';
    handles.pixelSize{3}.Enable = 'on';
  else
    handles.nPlanes.String = num2str(1);  
    handles.nPlanes.Enable = 'off';
    handles.pixCoordText{3}.Enable = 'off';
    handles.pixelSize{3}.Enable = 'off';
  end
end

function saveCB(hObj,event)
%   if ~checkControls()
%     return
%   end
  updateMetadata();
  kitSaveJobset(jobset);
  uiresume(gcf);
end

function updateMetadata()

  md = jobset.metadata{iMov};
  
  md.nFrames = str2double(handles.nFrames.String);
  md.nChannels = handles.nChannels;
  md.nPlanes = str2double(handles.nPlanes.String)*md.nFrames*md.nChannels;
  md.frameSize(3) = str2double(handles.nPlanes.String);
  for iChan = 1:3
    md.wavelength(iChan) = str2double(handles.wavelength{iChan}.String)/1000;
  end
  % construct frameTime
  if strcmp(handles.timeLapse.String,'N/A')
    frameTime = zeros(md.frameSize(3),1);
  else
    timeLapse = str2double(handles.timeLapse.String);
    for iTime = 1:md.nFrames
      timeLapsed = (iTime-1)*timeLapse;
      for iFrame = 1:md.frameSize(3)
        frameTime(iFrame,iTime) = timeLapsed;
      end
    end
  end
  md.frameTime = frameTime;
  md.is3D = handles.is3D.Value;
  md.na = str2double(handles.numAper.String);
  for iCoord = 1:3
    md.pixelSize(iCoord) = str2double(handles.pixelSize{iCoord}.String)/1000;
  end
  md.validated = 1;
  jobset.metadata{iMov} = md;
  
  applyAll = handles.applyAll.Value;
  
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

end % kitValidateMetadata
