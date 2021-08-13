function job = kitExampleRefinementTest(spots,movie,job, ...
                                             channel,verbose,debug)
% Perform and test refinement of spots by gaussian mixture model fitting
%
%take output from kitExampleDetectionTest.m rather than redoing
%reading in
% Jonathan U Harrison 2019-02-07
%%%%%%%%%%


%%%%% can get output from metaphase example by rerunning previous test
if nargin <1
    [spots, movie, job] = kitExampleDetectionTest();
end
if nargin<4 || isempty(channel)
channel = 1;
end
if nargin<5 || isempty(verbose)
    verbose=0;
end
if nargin<6 || isempty(debug)
    debug = 0;
end
method = 'gaussian'; %could also test centroid fitting method
ndims = 3;

% set options for use later in MMF fitting
job.options.debug.showMmfFinal = debug;
job.options.debug.showMmfPvals = debug;
job.options.debug.mmfVerbose = verbose;
job.options.spotMode{1} = 'adaptive';
job.options.chrShift.result{1,1} = zeros(1,6);
job.options.coordSystemChannel = 1;
%job.options.deconvolvedDataCorrection = 0; %change psf for deconvolved data

%get data properties rather than setting manually for test
ds = kitMakeMakiDatastruct(job,channel);
job.dataStruct{channel}.dataProperties = ds.dataProperties;

%default optins for mmf fitting
mmf.clusterSeparation = 5; % in PSF sigmas. If too large will fit whole
                            % plate, if too small will not account for
                            % overlapping PSFs.
mmf.alphaF = [0.05 0.05 0.05 0.05]; % N vs N+1 F-test cutoff.
mmf.alphaA = [0.05 0.075 0.10 0.10]; % amplitude t-test cutoff.
mmf.alphaD = [0.01 0.01 0.01 0.01]; % distance t-test cutoff.
mmf.mmfTol = 1e-5; % accuracy to which MMF fits Gaussians.
mmf.oneBigCluster = 0; % fit all spots together as one big cluster
mmf.maxMmfTime = 300; % Maximum per-frame time (sec) to attempt mixture model fit
                      % before giving up.  Use zero to disable.
mmf.addSpots = 0; % Try fitting multiple Gaussians to spots to identify overlaps
job.options.mmf = mmf;

% job.dataStruct{channel}.dataProperties.FILTERPRM(1:3) = job.dataStruct{channel}.dataProperties.FILTERPRM(1:3)/2;
% job.dataStruct{channel}.dataProperties.FT_SIGMA = job.dataStruct{channel}.dataProperties.FT_SIGMA/2;
%compute filters to use for background etc
filters = createFilters(ndims,job.dataStruct{channel}.dataProperties);

[imageSizeX,imageSizeY,imageSizeZ,nFrames] = size(movie);
localMaxima = repmat(struct('cands',[]),nFrames,1);
nSpots = zeros(nFrames,1);
for i=1:nFrames
    nSpots(i) = size(spots{i},1);
    
    % Round spots to nearest pixel and limit to image bounds.
    if nSpots(i) > 1
        spots{i} = bsxfun(@min,bsxfun(@max,round(spots{i}),1),[imageSizeX,imageSizeY,imageSizeZ]);
    end
    
    % Store the cands of the current image
    % TODO this is computed in both spot detectors, just return it.
    img = movie(:,:,:,i);
    if verLessThan('images','9.2')
        background = fastGauss3D(img,filters.backgroundP(1:3)/2,filters.backgroundP(4:6));
    else
        background = imgaussfilt3(img,filters.backgroundP(1:3)/2,'FilterSize',filters.backgroundP(4:6));
    end
    localMaxima(i).cands = spots{i};
    if nSpots(i) > 1
        spots1D = sub2ind(size(img),spots{i}(:,1),spots{i}(:,2),spots{i}(:,3));
    else
        spots1D = [];
    end
    localMaxima(i).candsAmp = img(spots1D);
    localMaxima(i).candsBg = background(spots1D);
    
    % Visualize candidates.
    if verbose ~= 0
        showSpots(img,spots{i});
        title(['Local maxima cands n=' num2str(size(spots{i},1))]);
        drawnow;
    end
end
kitLog('Average particles per frame: %.1f +/- %.1f',mean(nSpots),std(nSpots));

% Refine spot candidates.
switch method
    case 'centroid'
        job = kitCentroid(job,movie,localMaxima,channel);
    case 'gaussian'
        job = kitMixtureModel(job,movie,localMaxima,channel);
    case 'norefine'
        % No refinement. Copy localMaxima to initCoords.
        initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0, ...
            'amp',[],'bg',[]);
        initCoord(1).localMaxima = localMaxima;
        for i=1:nFrames
            initCoord(i).nSpots = size(localMaxima(i).cands,1);
            initCoord(i).allCoordPix = [localMaxima(i).cands(:,[2 1 3]) ...
                0.25*ones(initCoord(i).nSpots,3)];
            initCoord(i).allCoord = bsxfun(@times, initCoord(i).allCoordPix,...
                repmat(job.metadata.pixelSize,[1 2]));
            initCoord(i).amp = [localMaxima(i).candsAmp zeros(initCoord(i).nSpots,1)];
        end
        % Store data.
        job.dataStruct{channel}.initCoord = initCoord;
        job.dataStruct{channel}.failed = 0;
    otherwise
        error(['Unknown coordinate finding mode: ' method]);
end

nSpots = zeros(nFrames,1);
for i=1:nFrames
    nSpots(i) = job.dataStruct{channel}.initCoord(i).nSpots;
end    
fprintf('\nAverage number of refined spots in now: %.1f +/- %.1f\n',...
    mean(nSpots), ...
    std(nSpots));

%% Test1: check if reasonable number of spots
realisticNumSpots = 50;
assert(min(nSpots)>0.25*realisticNumSpots,'Should have at least a minimum number of spots')
assert(max(nSpots)<4*realisticNumSpots,'Should have no more than a max number of spots')
assert(size(nSpots,1)==nFrames,'Should have spots for each frame in movie')
assert(nSpots(1)>0, ...
    'Expect to find at least some kinetochores in the first frame')



