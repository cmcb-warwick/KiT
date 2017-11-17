function job=kitDefaultOptions()
% KITDEFAULTOPTIONS Populate default options struct
%
%    JOB = KITDEFAULTOPTIONS() Create job struct populated with default options.
%
% Copyright (c) 2013 Jonathan Armond

maxWavelengths = 3;

job.kit = 1;
job.version = kitVersion();
job.jobsetVersion = kitVersion(2);
job.software = ['KiT ' job.version];
job.matlabVersion = version;

job.movieDirectory = [];
job.movieFiles = [];

opts.coordMode{1} = 'centroid'; % Possible values: centroid, gaussian, none
opts.coordMode{2} = 'none';
opts.coordMode{3} = 'none';
opts.coordSystem = 'plate'; % Possible values: plate, image
opts.coordSystemChannel = 1;
opts.spotMode{1} = 'histcut'; % Possible values: histcut, wavelet, none
opts.spotMode{2} = 'none';
opts.spotMode{3} = 'none';
opts.cropAsked = 0;
opts.disableSave = 0; % Don't save job at each step. For debugging.

% Debug options.
debug.showMmfClusters = 0; % visualize clustering from overlapping PSFs, -1
                           % to pause, -2 to debug.
debug.showMmfCands = 0; % visualize centroid candidates, -1 to pause
debug.showMmfFinal = 0; % visualize final MMF spots, -1 to pause
debug.showMmfPvals = 0; % histogram of mixture-model p-values
debug.mmfVerbose = 0; % display detailed progress of mixture-model fitting.
debug.gapClosing = 0; % histogram of track gap lengths
debug.groupSisters = 0; % 1 - plot 4 frames with sisters assigment
                        % 2 - plot all tracks
debug.showIntensityMasks = 0;
debug.showPlaneFit = 0; % 1 to show plane fits, 2 to show each frame
debug.showCentroidFinal = 0; % visualize centroid final spots.
debug.showWavelet = 0; % 1 to show wavelet algorithm stages, 2 to save images.
debug.showWaveletAdapt = 0; % 1 to show adaptation of wavelet threshold.
debug.showAdaptive = 0; % 1 to show adaptive thresholding algorithm verbosely.
debug.asserts = 0; % check things that shouldn't go wrong
opts.debug = debug;


% Maki options.
opts.minSpotsPerFrame = 20;
opts.maxSpotsPerFrame = 1000;
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

% Gaussian mixture-model spot finding options.
opts.clusterSeparation = 5; % in PSF sigmas. If too large will fit whole
                            % plate, if too small we not account for
                            % overlapping PSFs.
opts.alphaF = [0.05 0.05  0.05]; % N vs N+1 F-test cutoff.
opts.alphaA = [0.05 0.075 0.10]; % amplitude t-test cutoff.
opts.alphaD = [0.01 0.01  0.01]; % distance t-test cutoff.
opts.oneBigCluster = 0; % fit all spots together as one big cluster
opts.maxMmfTime = 300; % Maximum per-frame time (sec) to attempt mixture model fit
                      % before giving up.  Use zero to disable.
opts.mmfAddSpots = 0; % Try fitting multiple Gaussians to spots to identify overlaps

% Image moment coordinate system.
opts.momentPrctile = 99; % percentile to threshold image at before computing
                         % moments

% Local spot intensity.
opts.maskRadius = 0.3; % um
opts.photobleachCorrect = 1;
opts.attachmentCorrect = 1;
opts.dirAssignMode = 'voting'; % directional (AP or P) assignment
opts.dirAssignExpWeight = 1; % exponentially weight displacements
opts.maskShape = 'circle';
opts.maskConeAngle = 10; % degrees, only used with maskShape=='cone'.
opts.poleShift = 0; % um
opts.otherSpotSearchRadius = 0;
opts.gaussFilterSpots = 0;

% Adaptive threshold spot detection.
opts.adaptiveLambda = 10; % Regularization for increasing number of spots.

% Deconvolution.
opts.deconvolve = 0;

% Wavelet multiscale product spot detector options.
opts.waveletLevelThresh = 2; % threshold scale for local MAD thresholding
opts.waveletLevelAdapt = 1; % use adaptive setting for above.
opts.waveletNumLevels = 3;  % number of wavelet levels
opts.waveletLocalMAD = 0; % locally estimated MAD
opts.waveletBackSub = 0;  % background subtraction
opts.waveletMinLevel = 1; % discard wavelet levels below this

% Manual spot detection options.
manualDetect.warningDist = 0.2; % in um, distance between spots below which to warn user
manualDetect.numFrames = 4; % minimum number of frames to manually select spots
manualDetect.gapMethod = 'framewise'; % only current option, frame-wise coming
opts.manualDetect = manualDetect;

job.options = opts;
