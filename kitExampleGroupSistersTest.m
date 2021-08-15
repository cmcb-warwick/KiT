function job = kitExampleGroupSistersTest(job,allowMulti,verbose,movie,theta,channel)
% Perform and test tracking
%
% take output from kitExampleDetectionTest.m and
% kitExampleRefinementTest.m and kitExampleFitTest.m and
% kitExampleGroupSistersTest.m
% rather than redoing
%
% movie input optional for visualisation
% Jonathan U Harrison 2019-02-12
%%%%%%%%%%

%%%%% can get output from metaphase example by rerunning previous tests
if nargin <1
    [spots, movie, job] = kitExampleDetectionTest();
    job = kitExampleRefinementTest(spots,movie,job);
    job = kitExamplePlaneFitTest(job);
    job = kitExampleTrackTest(job);
end

if nargin<2
    allowMulti=0;
end

if nargin<3 || isempty(verbose)
    verbose = 1;
end

if nargin<4
    movie = []; %then ignore the movie argument
end

if nargin<5
    %set default options for grouping sisters
    opts.useAlignment = 1;
    opts.maxAngle = 30;
    opts.maxDist = 1.5;
    opts.minOverlap = 10;
    opts.useAnaphase=0;
    opts.robust=0;
else %use input, can set up as optimization problem
    opts.useAlignment = 1;
    opts.maxAngle = theta(1);
    opts.maxDist = theta(2);
    opts.minOverlap = theta(3);
    opts.useAnaphase=0;
    opts.robust=0;
end
if nargin < 6
 channel = 1;
end

job.dataStruct{channel}.dataProperties.groupSisters=opts;

job.dataStruct{channel} = kitGroupSisters(job.dataStruct{channel},...
    verbose,allowMulti);

%how many sisters does this give us on average through the movie?
nFrames = job.metadata.nFrames;
nSisters = size(job.dataStruct{channel}.sisterList,1);
doesSisterPairExist = zeros(nFrames,nSisters);
if ~isempty(job.dataStruct{channel}.sisterList(1).coords1)
    for j = 1:nSisters
        %put ones where a paired sister exists
        doesSisterPairExist(:,j)= ~isnan(...
            job.dataStruct{channel}.sisterList(j).coords1(:,1));
    end
end
%multiply by 2 since we count the pair twice
nSistersMovieAverage = 2*mean(sum(doesSisterPairExist,2));
fprintf('On average through the movie we find %f sisters\n',...
    nSistersMovieAverage);

%can we plot which sisters we have matched up
if verbose && ~isempty(movie)
    frameToShow = min(50,nFrames); %pick a frame to plot sisters on
    kitVisualiseTrackedAndGroupedSpots(job,movie,channel,frameToShow);
elseif verbose
    %load movie if input movie is empty. eg:
    [~, movie, ~] = kitExampleDetectionTest(job);
    frameToShow = min(50,nFrames); %pick a frame to plot sisters on
    kitVisualiseTrackedAndGroupedSpots(job,movie,channel,frameToShow);
end
