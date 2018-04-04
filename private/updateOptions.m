function jobset = updateOptions(jobset,verbose)
% Update any older versions of a jobset's options
%

opts = jobset.options;
newOpts = kitDefaultOptions;
newOpts = newOpts.options;

if nargin<2
  verbose = 0;
end

if verbose
  kitLog('Current jobset version: %i',jobset.jobsetVersion)
  kitLog('Updating jobset to version %i.',kitVersion(2))
end
jobset.jobsetVersion = kitVersion(2);

newOpts.jobProcess = opts.jobProcess;
newOpts.spotMode = opts.spotMode;
newOpts.coordMode = opts.coordMode;
newOpts.coordSystem = opts.coordSystem;
newOpts.coordSystemChannel = opts.coordSystemChannel;
if isfield(opts,'cropAsked');
  newOpts.cropAsked = opts.cropAsked;
end

% Chromatic shift information.
if iscell(opts.chrShift)
  chrShift = newOpts.chrShift;
  chrShift.result = opts.chrShift;
  if verbose; kitLog('chrShift sub-structure generated'); end
else
  chrShift = opts.chrShift;
  % include coordinateAdjustments and regionFilter, if necessary
  if ~isfield(chrShift,'coordinateAdjustments')
    chrShift.coordinateAdjustments = repmat({zeros(1,3)},3,3);
    if verbose; kitLog('Included chrShift.coordinateAdjustments'); end
  end
  if ~isfield(chrShift,'regionFilter')
    chrShift.regionFilter = 1;
    if verbose; kitLog('Included chrShift.regionFilter'); end
  end
  % convert nnDist to neighbour, and amplitude to intensity
  if isfield(chrShift,'nnDistFilter')
    chrShift.neighbourFilter = chrShift.nnDistFilter;
    chrShift = rmfield(chrShift,'nnDistFilter');
    if verbose; kitLog('Converted chrShift.nnDistFilter to .neighbourFilter'); end
  end
  if isfield(chrShift,'amplitudeFilter')
    chrShift.intensityFilter = chrShift.amplitudeFilter;
    chrShift = rmfield(chrShift,'amplitudeFilter');
    if verbose; kitLog('Converted chrShift.amplitudeFilter to .intensityFilter'); end
  end
  % remove maskRadius, maskShape and interphase, if necessary
  if isfield(chrShift,'maskRadius')
    chrShift = rmfield(chrShift,'maskRadius');
    if verbose; kitLog('Removed chrShift.maskRadius'); end
  end
  if isfield(chrShift,'maskShape')
    chrShift = rmfield(chrShift,'maskShape');
    if verbose; kitLog('Removed chrShift.maskShape'); end
  end
  if isfield(chrShift,'interphase')
    chrShift = rmfield(chrShift,'interphase');
    if verbose; kitLog('Removed chrShift.interphase'); end
  end
end
newOpts.chrShift = chrShift;

% Debug options.
if isfield(opts,'debug')
  debug = opts.debug;
  if ~isfield(debug,'disableSave')
    debug.disableSave = 0;
  end
else
  debug = newOpts.debug;
end
newOpts.debug = debug;

% Maki options.
newOpts.minSpotsPerFrame = opts.minSpotsPerFrame;
if isfield(opts,'maxSpotsPerFrame')
  newOpts.maxSpotsPerFrame = opts.maxSpotsPerFrame;
else
  newOpts.maxSpotsPerFrame = 1000;
  if verbose; kitLog('Included maxSpotsPerFrame'); end
end
newOpts.betterBackground = opts.betterBackground;
newOpts.fitCloud = opts.fitCloud;
newOpts.robustStats = opts.robustStats;
newOpts.maxCloseGap = opts.maxCloseGap;
newOpts.autoRadiidt = opts.autoRadiidt;
newOpts.minSearchRadius = opts.minSearchRadius;
newOpts.maxSearchRadius = opts.maxSearchRadius;
newOpts.useSisterAlignment = opts.useSisterAlignment;
newOpts.maxSisterAlignmentAngle = opts.maxSisterAlignmentAngle;
newOpts.maxSisterSeparation = opts.maxSisterSeparation;
newOpts.minSisterTrackOverlap = opts.minSisterTrackOverlap;
newOpts.oppositeAnaphaseDir = opts.oppositeAnaphaseDir;

% Direction options.
if isfield(opts,'direction')
  direction = opts.direction;
else
  direction = newOpts.direction;
  if verbose; kitLog('Direction sub-structure generated'); end
  direction.assignMode = opts.dirAssignMode;
  direction.assignExpWeight = opts.dirAssignExpWeight;
end
newOpts.direction = direction;

% Gaussian mixture-model spot finding options.
if isfield(opts,'mmf')
  mmf = opts.mmf;
else
  mmf = newOpts.mmf;
  if verbose; kitLog('MMF sub-structure generated'); end
  mmf.clusterSeparation = opts.clusterSeparation;
  if length(opts.alphaF)==1
    mmf.alphaF = [opts.alphaF 0.05 0.05];
    if verbose; kitLog('Included channel-wise alphaF'); end
  else
    mmf.alphaF = opts.alphaF;
  end
  if length(opts.alphaA)==1
    mmf.alphaA = [opts.alphaA 0.075 0.10];
    if verbose; kitLog('Included channel-wise alphaA'); end
  else
    mmf.alphaA = opts.alphaA;
  end
  if length(opts.alphaD)==1
    mmf.alphaD = [opts.alphaD 0.01 0.01];
    if verbose; kitLog('Included channel-wise alphaD'); end
  else
    mmf.alphaD = opts.alphaD;
  end
  mmf.mmfTol = opts.mmfTol;
  mmf.oneBigCluster = opts.oneBigCluster;
  mmf.maxMmfTime = opts.maxMmfTime;
  mmf.addSpots = opts.mmfAddSpots;
end
newOpts.mmf = mmf;

% Image moment coordinate system.
newOpts.momentPrctile = opts.momentPrctile;

% Neighbouring spot localisation.
if isfield(opts,'neighbourSpots')
  neighbourSpots = opts.neighbourSpots;
  if ~isfield(neighbourSpots,'channelOrientation');
    neighbourSpots.channelOrientation = 1:3;
    if verbose; kitLog('Included neighbourSpots.channelOrientation'); end
  end
else
  neighbourSpots = newOpts.neighbourSpots;
  if verbose; kitLog('NeighbourSpots sub-structure generated'); end
end
newOpts.neighbourSpots = neighbourSpots;

% Local spot intensity.
if isfield(opts,'intensity')
  intensity = opts.intensity;
  if isfield(intensity,'do')
    intensity.execute = intensity.do;
    intensity = rmfield(intensity,'do');
    if verbose; kitLog('Converted intensity.do to .execute'); end
  end
else
  intensity = newOpts.intensity;
  if verbose; kitLog('Intensity sub-structure generated'); end
  intensity.maskRadius = opts.maskRadius;
  intensity.photobleachCorrect = opts.photobleachCorrect;
  intensity.attachmentCorrect = opts.attachmentCorrect;
  intensity.maskShape = opts.maskShape;
  intensity.maskConeAngle = opts.maskConeAngle;
  intensity.poleShift = opts.poleShift;
  intensity.otherSpotSearchRadius = opts.otherSpotSearchRadius;
  intensity.gaussFilterSpots = opts.gaussFilterSpots;
end
newOpts.intensity = intensity;

% Adaptive threshold spot detection.
newOpts.adaptiveLambda = opts.adaptiveLambda;

% Wavelet multiscale product spot detector options.
if isfield(opts,'wavelet')
  wavelet = opts.wavelet;
else
  wavelet = newOpts.wavelet;
  if verbose; kitLog('Wavelet sub-structure generated'); end
  wavelet.levelThresh = opts.waveletLevelThresh;
  wavelet.levelAdapt = opts.waveletLevelAdapt;
  wavelet.numLevels = opts.waveletNumLevels;
  wavelet.localMAD = opts.waveletLocalMAD;
  wavelet.backSub = opts.waveletBackSub;
  wavelet.minLevel = opts.waveletMinLevel;
end
newOpts.wavelet = wavelet;

% Manual spot detection options.
if isfield(opts,'manualDetect')
  manualDetect = opts.manualDetect;
else
  manualDetect = newOpts.manualDetect;
  if verbose; kitLog('ManualDetect sub-structure generated'); end
end
newOpts.manualDetect = manualDetect;

jobset.options = newOpts;

