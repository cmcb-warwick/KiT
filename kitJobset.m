function jobset=kitJobset(jobset)
% KITJOBSET Create or upgrade a jobset.

if nargin<1
  jobset = kitDefaultOptions();
  return
end

% Upgrade to use ROIs.
if jobset.jobsetVersion < 7
  jobset = updVer7(jobset);
end

% Modularisation of options.
if jobset.jobsetVersion < 8
  % Update options.
  jobset.options = updVer8(jobset.options);
end

% Upgrade to using 4 channels.
if jobset.jobsetVersion < 9
  jobset.options = updVer9(jobset.options);
end

% Copy any missing fields.
defJob = kitDefaultOptions();
jobset = structCopyMissingFields(jobset,defJob);

% Update jobset version.
jobset.jobsetVersion = kitVersion(2);

end

function js = updVer7(js)
  % Upgrade to use ROIs.
  for i=1:length(js.ROI)
    js.ROI(i).movie = js.movieFiles{i};
    js.ROI(i).crop = js.crop{i};
    js.ROI(i).cropSize = js.cropSize{i};
  end
  js = rmfield(js,{'crop','cropSize'});
  js.movieFiles = unique(js.movieFiles);
  
end

function opts = updVer8(opts)
  % Modularisation of options.
  
  % Get default options.
  defOpts = kitDefaultOptions;
  defOpts = defOpts.options;

  % Chromatic shift information.
  if iscell(opts.chrShift)
    cS = defOpts.chrShift;
    cS.result = opts.chrShift;
  else
    cS = opts.chrShift;
    % convert nnDist to neighbour, and amplitude to intensity
    if isfield(cS,'nnDistFilter')
      cS.neighbourFilter = cS.nnDistFilter;
      cS = rmfield(cS,'nnDistFilter');
    end
    if isfield(cS,'amplitudeFilter')
      cS.intensityFilter = cS.amplitudeFilter;
      cS = rmfield(cS,'amplitudeFilter');
    end
    % remove maskRadius, maskShape and interphase, if necessary
    if isfield(cS,'maskRadius')
      cS = rmfield(cS,'maskRadius');
    end
    if isfield(cS,'maskShape')
      cS = rmfield(cS,'maskShape');
    end
    if isfield(cS,'interphase')
      cS = rmfield(cS,'interphase');
    end
  end
  opts.chrShift = cS;

  % Debug option.
  if isfield(opts,'disableSave')
    opts.debug.disableSave = 0;
    opts = rmfield(opts,'disableSave');
  end

  % Direction options.
  if ~isfield(opts,'direction')
    dir = defOpts.direction;
    dir.assignMode = opts.dirAssignMode;
    dir.assignExpWeight = opts.dirAssignExpWeight;
    opts.direction = dir;
    opts = rmfield(opts,{'dirAssignMode','dirAssignExpWeight'});
  end

  % Gaussian mixture-model spot finding options.
  if ~isfield(opts,'mmf')
    mmf = defOpts.mmf;
    if length(opts.alphaA)==1
      mmf.alphaA = [opts.alphaA 0.075 0.10];
    else
      mmf.alphaA = opts.alphaA;
    end
    if length(opts.alphaD)==1
      mmf.alphaD = [opts.alphaD 0.01 0.01];
    else
      mmf.alphaD = opts.alphaD;
    end
    if length(opts.alphaF)==1
      mmf.alphaF = [opts.alphaF 0.05 0.05];
    else
      mmf.alphaF = opts.alphaF;
    end
    mmf.clusterSeparation = opts.clusterSeparation;
    mmf.mmfTol = opts.mmfTol;
    mmf.oneBigCluster = opts.oneBigCluster;
    mmf.maxMmfTime = opts.maxMmfTime;
    mmf.addSpots = opts.mmfAddSpots;
    opts = rmfield(opts,{'mmfTol','oneBigCluster','maxMmfTime','mmfAddSpots'});
  end
  opts.mmf = mmf;

  % Neighbouring spot localisation.
  if ~isfield(opts.neighbourSpots,'channelOrientation')
    opts.neighbourSpots.channelOrientation = 1:3;
  end

  % Local spot intensity.
  if isfield(opts,'intensity')
    if isfield(opts.intensity,'do')
      opts.intensity.execute = opts.intensity.do;
      opts.intensity = rmfield(opts.intensity,'do');
    end
  else
    intensity = defOpts.intensity;
    intensity.maskRadius = opts.maskRadius;
    intensity.photobleachCorrect = opts.photobleachCorrect;
    intensity.attachmentCorrect = opts.attachmentCorrect;
    intensity.maskShape = opts.maskShape;
    intensity.maskConeAngle = opts.maskConeAngle;
    intensity.poleShift = opts.poleShift;
    intensity.otherSpotSearchRadius = opts.otherSpotSearchRadius;
    intensity.gaussFilterSpots = opts.gaussFilterSpots;
    opts = rmfield(opts,{'maskRadius','photobleachCorrect','attachmentCorrect',...
        'maskShape','maskConeAngle','poleShift','otherSpotSearchRadius',...
        'gaussFilterSpots'});
    opts.intensity = intensity;
  end

  % Wavelet multiscale product spot detector options.
  if ~isfield(opts,'wavelet')
    wavelet = defOpts.wavelet;
    wavelet.levelThresh = opts.waveletLevelThresh;
    wavelet.levelAdapt = opts.waveletLevelAdapt;
    wavelet.numLevels = opts.waveletNumLevels;
    wavelet.localMAD = opts.waveletLocalMAD;
    wavelet.backSub = opts.waveletBackSub;
    wavelet.minLevel = opts.waveletMinLevel;
    opts = rmfield(opts,{'waveletLevelThresh','waveletLevelAdapt','waveletNumLevels',...
        'waveletLocalMAD','waveletBackSub','waveletMinLevel'});
  end

end

function opts = updVer9(opts)
%   % Upgrade to 4 channel.  
% 
%   % Detection and refinement modes.
%   opts.spotMode{4} = 'none';
%   opts.coordMode{4} = 'none';
%   
%   % Chromatic shift.
%   cSresult = repmat({zeros(1,6)},4,4);
%   cSresult(1:3,1:3) = opts.chrShift.result;
%   opts.chrShift.result = cSresult;
%   cSjobset = repmat({[]},4,4);
%   cSjobset(1:3,1:3) = opts.chrShift.jobset;
%   opts.chrShift.jobset = cSjobset;
%   for i=1:4; for j=1:4; chanOrder{i,j} = [i j]; end; end
%   chanOrder(1:3,1:3) = opts.chrShift.chanOrder;
%   opts.chrShift.chanOrder = chanOrder;
%   coordAdjs = repmat({zeros(1,3)},4,4);
%   coordAdjs(1:3,1:3) = opts.chrShift.coordinateAdjustments;
%   opts.chrShift.coordinateAdjustments = coordAdjs;
%   
%   % MMF.
%   opts.mmf.alphaA(4) = 0.10;
%   opts.mmf.alphaD(4) = 0.01;
%   opts.mmf.alphaF(4) = 0.05;
%   
%   % Neighbour spot detection.
%   opts.neighbourSpots.channelOrientation(4) = 4;
%   opts.neighbourSpots.timePoints{4} = [];
%   opts.neighbourSpots.zSlices{4} = [];
%   
%   % Intensity measurement.
%   opts.intensity.execute(4) = 0;
%   
end