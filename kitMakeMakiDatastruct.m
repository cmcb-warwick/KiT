function dataStruct = kitMakeMakiDatastruct(job, channel)
%KITMAKEMAKIDATASTRUCT Generates maki compatible datastruct from kit job
%
% SYNOPSIS: dataStruct = kitMakeMakiDatastruct(job)
%
% INPUT job: Struct containing tracking job setup options.
%            Requires at least the following fields:
%
%           .
%
% OUTPUT dataStruct: Maki compatible struct containing the following fields:
%
%       .rawMovieName - may contain the actual movie
%       .rawMoviePath
%       .dataProperties - struct containing
%
% Copyright (c) 2012 Jonathan W. Armond

options = job.options;

% Default values for non-Kit fields.
dataProperties.LENSID=12003;
dataProperties.timeLapse=1;
dataProperties.expTime=NaN;
dataProperties.NDfilter=NaN;
dataProperties.cellCycle=NaN;
dataProperties.strains=NaN;
dataProperties.drugs=NaN;
dataProperties.temperature={'Nan'};
dataProperties.crop=[];
dataProperties.F_TEST_PROB=0.9990;
dataProperties.IDopt= [];
dataProperties.PATCHSIZE=7;
dataProperties.CH_MAXNUMINTERV=1000;
dataProperties.OVERLPSIZE=[15 15 15];
dataProperties.sigmaCorrection=[1.5 1.5];
dataProperties.split=[];
dataProperties.MAXSPOTS=500;
dataProperties.T_TEST_PROB=0.0500;
dataProperties.maxSize = 100 * 2^20; % 100 Mb
dataProperties.amplitudeCutoff = 0; % undefined
dataProperties.fitNPlusOne = 1; % super-resolution fitting in detector
dataProperties.waveIdx = channel; % current wavelength
dataProperties.movieType = 'sorger'; % also: 'sedat', 'misteli', 'synth'
dataProperties.name = '';

% linker properties
dataProperties.linker_relativeMaxDistance = -1; % don't use
dataProperties.linker_absoluteMaxDistance=-1; % don't use
dataProperties.linker_relAmpWeight=1/1.5; % weighs distance more
dataProperties.linker_useCOM = 1; % use center of mass to correct
dataProperties.linker_fuseRatio = 1.5; % fuse if less than 1.5 RL separated

% detector properties
dataProperties.detector_spotfind = 1; %1: standard 2: mammalian

% tracking settings
tracksParam.rotate = 1;
tracksParam.timeWindow = options.maxCloseGap;
tracksParam.minRadius = options.minSearchRadius;
tracksParam.maxRadius = options.maxSearchRadius;

% gap closing parameters.
gapCloseParam.timeWindow = tracksParam.timeWindow + 1;
gapCloseParam.mergeSplit = 0;
gapCloseParam.minTrackLen = 2;
gapCloseParam.diagnostics = options.debug.gapClosing;

%assign cost matrix parameters for linking spots between consecutive
%frames
costMatrices(1).funcName = 'trackCostMatLink';
parameters.linearMotion = 0;
parameters.minSearchRadius = tracksParam.minRadius;
parameters.maxSearchRadius = tracksParam.maxRadius;
parameters.brownStdMult = 3.5;
parameters.useLocalDensity = 1;
parameters.nnWindow = gapCloseParam.timeWindow;
costMatrices(1).parameters = parameters;
clear parameters

%assign cost matrix parameters for closing gaps and (in principle)
%merging and splitting
costMatrices(2).funcName = 'trackCostMatCloseGaps';
parameters.linearMotion = 0;
parameters.minSearchRadius = tracksParam.minRadius;
parameters.maxSearchRadius = tracksParam.maxRadius;
parameters.brownStdMult = 3.5*ones(gapCloseParam.timeWindow,1);
% parameters.timeReachConfB = min(2,gapCloseParam.timeWindow);
parameters.timeReachConfB = min(1,gapCloseParam.timeWindow);
parameters.lenForClassify = 10;
parameters.ampRatioLimit = [0.65 4];
parameters.useLocalDensity = 1;
parameters.nnWindow = gapCloseParam.timeWindow;
parameters.linStdMult = 3.5*ones(gapCloseParam.timeWindow,1);
parameters.timeReachConfL = 1;
parameters.maxAngleVV = 45;
costMatrices(2).parameters = parameters;
clear parameters

%assign Kalman filter function names
kalmanFunctions.reserveMem = 'kalmanResMemLM';
kalmanFunctions.initialize = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%save tracking parameters in dataStruct
tracksParam.gapCloseParam = gapCloseParam;
tracksParam.costMatrices = costMatrices;
tracksParam.kalmanFunctions = kalmanFunctions;
dataProperties.tracksParam = tracksParam;

% Group sisters.
groupSisters.useAlignment = options.useSisterAlignment;
groupSisters.maxAngle = options.maxSisterAlignmentAngle;
groupSisters.maxDist = options.maxSisterSeparation;
groupSisters.minOverlap = options.minSisterTrackOverlap;
groupSisters.useAnaphase = options.oppositeAnaphaseDir;
groupSisters.robust = options.robustStats;
dataProperties.groupSisters = groupSisters;


% Translate metadata.
md = job.metadata;

dataProperties.PIXELSIZE_XY = md.pixelSize(1);
if md.pixelSize(1) ~= md.pixelSize(2)
    error('unequal pixelsize x/y!')
end
dataProperties.PIXELSIZE_Z = md.pixelSize(3);
dataProperties.NA = md.na;
if isempty(job.crop)
    dataProperties.movieSize(1:2) = md.frameSize(1:2);
else
    dataProperties.movieSize(1:2) = job.cropSize(1:2);
end
dataProperties.movieSize(3) = md.frameSize(3);
dataProperties.movieSize(4) = md.nFrames;
dataProperties.WVL = md.wavelength;
% if there are multiple wavelengths: Make it a numZ x numT x
% numW array
dataProperties.frameTime = repmat(md.frameTime, [1, 1, length(md.wavelength)]);

% calculate filterparms
[FT_XY, FT_Z] = calcFilterParms(...
  dataProperties.WVL(dataProperties.waveIdx),dataProperties.NA,[],'gauss',...
  dataProperties.sigmaCorrection,...
  [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z]);
patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];
dataProperties.FT_SIGMA = [FT_XY,FT_XY,FT_Z];


% Assign to dataStruct.
dataStruct.dataProperties = dataProperties;

