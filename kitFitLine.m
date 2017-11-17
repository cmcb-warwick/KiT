function dataStruct=kitFitLine(dataStruct)

if isfield(dataStruct,'planeFit')
    useInputPlane = 1;
else
    useInputPlane = 0;
end

% minimal number of consecutive frames in a movie that have stable enough
% eigenvectors for plane rotation estimation
minConsecFrames = 5; 

% get coordinates
initCoord = dataStruct.initCoord;
nTimePoints = length(initCoord);
nSpots = cat(1,initCoord.nSpots);

% setup planeFit structures.
if useInputPlane
    planeFit = dataStruct.planeFit;
else
    planeFit(1:nTimePoints) = struct('plane',[],'planeCoord',[],'planeVectorClassifier', 0, ...
        'planeVectors',[],'planeOrigin',[],'eigenVectors',[],'eigenValues',[],...
        'inlierIdx',[],'unalignedIdx',[],'laggingIdx',[],'phase','e',...
        'distParms',[],'deltaP',[],'deltaAngle',[]);
end 

% initialize lists of frames with and without plane
framesNoPlane = 1:nTimePoints;
framesWiPlane = [];

% probDim 

% loop through timepoints. Get covariance of point cloud, and the
% corresponding eigenvalues. Label frames that have sufficient anisotropy.
eigenValues = zeros(nTimePoints,probDim);
eigenVectors = zeros(probDim,probDim,nTimePoints);  %vectors in cols
meanCoord = zeros(nTimePoints,probDim);
meanCoordFull = zeros(nTimePoints,3);

goodFrames1 = [];
potFrames = [];
for t=1:nTimePoints
