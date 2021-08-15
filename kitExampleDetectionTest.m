function [spots, movie, job] = kitExampleDetectionTest(job, channel,...
                            useAdaptive, adaptiveLambda, verbose, ignoreTests)
% Perform detection on metaphase example downloaded from the bitbucket repo
%
% Jonathan U Harrison 2019-02-06
%%%%%%%%%%
if nargin <1 || isempty(job)
    %use Metaphase Example movie
    %put together jobset struct so can load data etc
    job = kitDefaultOptions();
    job.movieDirectory = 'testdata';
    job.ROI.movie = 'MetaphaseExample.ome.tiff';
    job.ROI.crop = [1 1 109 134];
    job.ROI.cropSize = [109 134 25];
end

if ~isfield(job,'index')
    job.index=1;
end
%set some options for which method to use
if nargin<2 || isempty(channel)
    channel=1;
end
if nargin<3 || isempty(useAdaptive)
    useAdaptive=1;
end
if nargin<4 || isempty(adaptiveLambda) && useAdaptive
    adaptiveLambda = 1;
end
if nargin<5 || isempty(verbose)
    verbose=0;
end
if nargin<6 || isempty(ignoreTests)
    ignoreTests=0;
end
realisticNumSpots = job.options.realisticNumSpots;
globalBackground = job.options.globalBackground;

%create reader to read movie
%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(job,'ROI') && length(job.ROI)>1
    warning('More than one stack in job but using only first one\n');
    pathToMovie = fullfile(job.movieDirectory,job.ROI(1).movie);
else
    pathToMovie = fullfile(job.movieDirectory,job.ROI.movie);
end
if isfield(job,'metadata')
    if iscell(job.metadata)
        job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(pathToMovie,'valid',job.metadata);
else
    [job.metadata, reader] = kitOpenMovie(pathToMovie);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%load movie
movie = kitReadWholeMovie(reader,job.metadata,channel,job.ROI(1).crop,0,1);


if useAdaptive
    spots = adaptiveSpots(movie,adaptiveLambda,realisticNumSpots,globalBackground,1,0);
else
    %fix some parameters used in the histcut algorithm
    options.minSpotsPerFrame=80;
    options.maxSpotsPerFrame=110;
    dataProperties.FT_SIGMA = [0.1181 0.1181 0.4025];
    dataProperties.FILTERPRM= [0.1181 0.1181 0.4025 1 1 3];
    dataProperties.movieSize = size(movie);
    nFrames = size(movie,4);
    spots = cell(nFrames,1);
    for i=1:nFrames
        img = movie(:,:,:,i);
        spots{i} = histcutSpotsSimplified(img,options,dataProperties,verbose);
    end
end
%max project and show spots
if verbose
    t_frame = 1;
    figure;
    imshow(max(movie(:,:,:,t_frame),[],3));
    hold on;
    spots_max_project = max(spots{t_frame},[],3);
    scatter(spots_max_project(:,2),spots_max_project(:,1),'rx');
end

nSpotsDetected = zeros(job.metadata.nFrames,1);
for ii=1:job.metadata.nFrames
    nSpotsDetected(ii) = size(spots{ii},1);
end
figure; plot(1:job.metadata.nFrames,nSpotsDetected);
mean(nSpotsDetected)

%%%%%%%%%%%%%%%
% to write tests specify expectations of length of spots, and what spots
% should be
%%%%%%%%%%%%%%%%%%%
%% Test1: spots has the right dimensions
if ~ignoreTests
assert(length(spots)==size(movie,4),'Missing spots output for some frames')
assert(size(spots{1},2)==3,'Expecting a 3D movie to test on')
assert(size(spots{1},1)>0, ...
    'Expect to find at least some kinetochores in the first frame')

%% Test2: Found roughly the right amount of spots
nSpotsPerFrame = zeros(size(spots));
nFrames = size(movie,4);
for j = 1:nFrames
    nSpotsPerFrame(j) = size(spots{j},1);
end
assert(min(nSpotsPerFrame)>0.25*realisticNumSpots, ...
    sprintf('Expect a biologically reasonable number of spots per frame, found min of %d',min(nSpotsPerFrame)));
assert(max(nSpotsPerFrame)<4*realisticNumSpots, ...
    sprintf('Expect a biologically reasonable number of spots per frame, found max of %d',max(nSpotsPerFrame)));
%% Test3: does a spot seem like a genuine spot
assert(all(spots{1}(1,:)>0), 'Spot coords should be positive')
% assert(all(movie([spots{1}(1,:),1])>0), ...
%     'First spot in first frame should have nonzero amplitude')
% assert(all(movie([spots{1}(1,:),1])>median(median(median(median(movie))))), ...
%     'First spot in first frame should have amplitude grater than average')
end
