function job=kitDefaultOptions()
% KITDEFAULTOPTIONS Populate default options struct
%
%    JOB = KITDEFAULTOPTIONS() Create job struct populated with default options.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

job.kit = 1;
job.version = kitVersion();
job.jobsetVersion = kitVersion(2);
job.software = ['KiT ' job.version];
job.matlabVersion = version;

job.movieDirectory = [];
job.movieFiles = [];

opts.jobProcess = 'zandt';
opts.spotMode{1} = 'histcut'; % Possible values: histcut, wavelet, neighbour, manual, none
opts.spotMode{2} = 'none';
opts.spotMode{3} = 'none';
opts.coordMode{1} = 'gaussian'; % Possible values: centroid, gaussian, none
opts.coordMode{2} = 'none';
opts.coordMode{3} = 'none';
opts.coordSystem = 'plate'; % Possible values: plate, image, com
opts.coordSystemChannel = 1;
opts.cropAsked = 0;

%Chromatic shift information
chrShift.result = repmat({zeros(1,6)},3,3);
chrShift.jobset = repmat({[]},3,3);
for i=1:3; for j=1:3; chrShift.chanOrder{i,j} = [i j]; end; end
chrShift.coordinateAdjustments = repmat({zeros(1,3)},3,3);
chrShift.minSpots = 20;
chrShift.filtering = 0;
chrShift.intensityFilter = 0; % ratio of spot intensity with max intensity
chrShift.neighbourFilter = 0; % in um, min distance between spots
chrShift.regionFilter = 1; % number of regions to divide chromatic shift images
opts.chrShift = chrShift;

% Debug options.
debug.showMmfClusters = 0; % visualize clustering from overlapping PSFs, -1
                           % to pause, -2 to debug.
debug.showMmfCands = 0; % visualize centroid candidates, -1 to pause
debug.showMmfFinal = 0; % visualize final MMF spots, -1 to pause
debug.showMmfPvals = 0; % histogram of mixture-model p-values
debug.mmfVerbose = 0; % display detailed progress of mixture-model fitting.
debug.gapClosing = 0; % histogram of track gap lengths
debug.groupSisters = 0; % 1 - plot 4 frames with sisters assignment
                        % 2 - plot all tracks
debug.showIntensityMasks = 0;
debug.showPlaneFit = 0; % 1 to show plane fits, 2 to show each frame
debug.showCentroidFinal = 0; % visualize centroid final spots.
debug.showWavelet = 0; % 1 to show wavelet algorithm stages, 2 to save images.
debug.showWaveletAdapt = 0; % 1 to show adaptation of wavelet threshold.
debug.showAdaptive = 0; % 1 to show adaptive thresholding algorithm verbosely.
debug.asserts = 0; % check things that shouldn't go wrong
debug.disableSave = 0; % don't save job at each step
opts.debug = debug;


% Maki options.
opts.minSpotsPerFrame = 20;
opts.maxSpotsPerFrame = 150;
opts.betterBackground = 0; % 1 == mask signal before taking background.
opts.fitCloud = 0; % 1 == use max evector as normal axis.
opts.robustStats = 0;
opts.maxCloseGap = 4;
opts.autoRadiidt = 2;
opts.minSearchRadius = [0.01 0.1 0.01]; % [inliers, unaligned, lagging] um
opts.maxSearchRadius = [0.75 3 0.75]; % [inliers, unaligned, lagging] um
opts.useSisterAlignment = 1;
opts.maxSisterAlignmentAngle = 30; % degrees
opts.maxSisterSeparation = 1.5; % um
opts.minSisterTrackOverlap = 10; % percent of movie length to require overlap
opts.oppositeAnaphaseDir = 1; % Use assumption of opposition direction in
                                 % anaphase.

direction.assignMode = 'voting'; % directional (AP or P) assignment, can also be 'absolute'
direction.assignExpWeight = 1; % exponentially weight displacements
direction.minConsSteps = 3; % minimum consecutive timepoints in one direction
direction.switchBuffer = 1; % number of timepoints at start and end of run
                          % to be considered directional switches
opts.direction = direction;

% Gaussian mixture-model spot finding options.
mmf.clusterSeparation = 5; % in PSF sigmas. If too large will fit whole
                            % plate, if too small will not account for
                            % overlapping PSFs.
mmf.alphaF = [0.05 0.05 0.05]; % N vs N+1 F-test cutoff.
mmf.alphaA = [0.05 0.075 0.10]; % amplitude t-test cutoff.
mmf.alphaD = [0.01 0.01 0.01]; % distance t-test cutoff.
mmf.mmfTol = 1e-5; % accuracy to which MMF fits Gaussians.
mmf.oneBigCluster = 0; % fit all spots together as one big cluster
mmf.maxMmfTime = 300; % Maximum per-frame time (sec) to attempt mixture model fit
                      % before giving up.  Use zero to disable.
mmf.addSpots = 0; % Try fitting multiple Gaussians to spots to identify overlaps
opts.mmf = mmf;

% Image moment coordinate system.
opts.momentPrctile = 99; % percentile to threshold image at before computing
                         % moments

% Neighbouring spot localisation.
neighbourSpots.maskRadius = 0.3; % um
neighbourSpots.maskShape = 'semicirc';
neighbourSpots.maskConeAngle = 10; % degrees, only used with maskShape=='cone'.
neighbourSpots.channelOrientation = 1:3; % from inner to outer
neighbourSpots.timePoints = repmat({[]},1,3);
neighbourSpots.zSlices = repmat({[]},1,3);
opts.neighbourSpots = neighbourSpots;

% Local spot intensity.
intensity.execute = zeros(1,3); % whether or not to run intensity measurements per channel
intensity.maskRadius = 0.3; % um
intensity.photobleachCorrect = 1;
intensity.attachmentCorrect = 1;
intensity.maskShape = 'circle';
intensity.maskConeAngle = 10; % degrees, only used with maskShape=='cone'.
intensity.poleShift = 0; % um
intensity.otherSpotSearchRadius = 0;
intensity.gaussFilterSpots = 0;
opts.intensity = intensity;

% Adaptive threshold spot detection.
opts.adaptiveLambda = 10; % Regularization for increasing number of spots.

% Deconvolution.
opts.deconvolve = 0;

% Wavelet multiscale product spot detector options.
wavelet.levelThresh = 2; % threshold scale for local MAD thresholding
wavelet.levelAdapt = 1; % use adaptive setting for above.
wavelet.numLevels = 3;  % number of wavelet levels
wavelet.localMAD = 0; % locally estimated MAD
wavelet.backSub = 0;  % background subtraction
wavelet.minLevel = 1; % discard wavelet levels below this
opts.wavelet = wavelet;

% Manual spot detection options.
manualDetect.warningDist = 0.2; % in um, distance between spots below which to warn user
manualDetect.frameSpacing = 7; % number of frames between adjacent manual detection frame
manualDetect.gapMethod = 'framewise'; % can also be 'linear'
opts.manualDetect = manualDetect;

job.options = opts;
