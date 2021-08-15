function kitDiagnoseAndTest(job,useExistingAnalysis,applyTests,spots,movie,channel)
%kitDiagnoseAndTest(job,useExistingAnalysis,applyTests,spots,movie)
%Keep track of how many sisters detected at each step.
%This may not necessarily be the best way to structure this diagnostic
%since may stop if a test fails?
%
%job - kit Jobset to specify the movie and analysis
%useExistingAnalysis - bool, run full tracking or use existing results.
%If using existing results, will still need to rerun the initial detection,
%as don't think this is stored.
%
%Pass empty argument [] to use metaphase example tiff
%Jonathan U harrison 2019-02-15
%%%%%%%%%%%%
if nargin<3
    applyTests = 0;
end    
if nargin<2
    useExistingAnalysis=0;
end
if nargin<1
    job=[]; %use metaphase example tiff by default
end
if nargin<4
%detect
    [spots, movie, job] = kitExampleDetectionTest(job);
%%%%%%%%%%%
end
if nargin<6
    channel=1;
end
if useExistingAnalysis && ~isfield(job,'dataStruct')
    error(['Tried to use existing analysis but looks' ...
        'like there is not one for this job']);
end
%%%%%%%%%%%
if ~useExistingAnalysis
    failed=0;
    while ~failed
        %refine via mmf
        try
            job = kitExampleRefinementTest(spots,movie,job);
        catch ME
            warning('failure at refinement step: %s \n',ME.identifier);
            failed=1;
        end
        %%%%%%%%%%%
        %fit plane
        try
            job = kitExamplePlaneFitTest(job);
        catch
            warning('failure at plane fitting step');
            failed=2;
        end
        %%%%%%%%%%%
        %track
        try
            job = kitExampleTrackTest(job);
        catch
            warning('failure at tracking step');
            failed=3;
        end
        %%%%%%%%%%%
        %group sisters
        try
            job = kitExampleGroupSistersTest(job);
        catch ME
            warning('failure at sister grouping step: %s\n',ME.identifier);
            failed=4;
        end
    end
else 
    failed = job.dataStruct{channel}.failed;
end
%%%%%%%%%%%
nFrames = size(spots,1);
nStages = 5; %number of stages of algorithm to consider
nSpots = zeros(nFrames,nStages); %store num spots at each different stage
if isfield(job,'dataStruct')
    minOverlap = job.dataStruct{channel}.dataProperties.groupSisters.minOverlap; %now in number of frames
else
    minOverlap = 10;
end
try
    planeFitted = zeros(nFrames,1);
    for i = 1:nFrames
        nSpots(i,1)=size(spots{i},1);
        if failed~=1 %even if failed at later stage should still have this
            nSpots(i,2)=job.dataStruct{channel}.initCoord(i).nSpots;
            if failed>2 || ~failed
                planeFitted(i) = (~isempty(job.dataStruct{channel}.planeFit)&&...
                    ~isempty(job.dataStruct{channel}.planeFit(i).planeVectors));
            end
        end
    end
    %compute how many tracks in a given frame
    nTracks = size(job.dataStruct{channel}.tracks,1);
    doesTrackExist = zeros(nFrames,nTracks);
    isTrackLong = zeros(nFrames,nTracks);
    for j = 1:nTracks
        %put ones where a track exists
        startTime = job.dataStruct{channel}.tracks(j).seqOfEvents(1,1);
        endTime = job.dataStruct{channel}.tracks(j).seqOfEvents(2,1);
        doesTrackExist(startTime:endTime,j) = ones(endTime-startTime+1,1);
        isTrackLong(startTime:endTime,j) = (endTime - startTime > ...
            minOverlap)*ones(endTime-startTime+1,1);
    end
    nSpots(:,3) = sum(doesTrackExist,2);
    nSpots(:,4) = sum(isTrackLong,2);
    %compute how many sisters in a given frame
    nSisters = size(job.dataStruct{channel}.sisterList,1);
    doesSisterPairExist = zeros(nFrames,nSisters);
    for j = 1:nSisters
        %put ones where a paired sister exists
        doesSisterPairExist(:,j)= ~isnan(...
            job.dataStruct{channel}.sisterList(j).coords1(:,1));
    end
    %multiply by 2 since we count the pair twice
    nSpots(:,5) = 2*sum(doesSisterPairExist,2);
catch
    warning(['Definitely did not complete properly and ' ...
        'will plot 0 for some stages']);
end
%%%%%%%%%%
%plot it all
t = job.metadata.frameTime(1,:);
if nStages ==5
    coloursByStage = ...
        [228,26,28;
        55,126,184;
        77,175,74;
        152,78,163;
        255,127,0]/256;
%         [127,201,127;
%         190,174,212;
%         253,192,134;
%         255,255,153;
%         56,108,176]/256; %from http://colorbrewer2.org/
elseif nStages==4
    % for 4 colours:
    coloursByStage = ...
       [230,97,1;
        253,184,99;
        178,171,210;
        94,60,153]/256; %from http://colorbrewer2.org/
else 
    error('Do not have colours for that number of stages');
end
symbolByPlane = 'xo';
stagesLegend = {'detect','mmf','track',...
    sprintf('>%d frames',minOverlap),'group'};
planeFitted = logical(planeFitted); %else 0s or 1s not treated properly
figure;
hold all;
for i=1:nStages
    if (sum(planeFitted)>0)
        scatter(t(planeFitted),nSpots(planeFitted,i),100,...
            repmat(coloursByStage(i,:),sum(planeFitted),1), symbolByPlane(1));
    elseif sum(planeFitted)<nFrames
        warning('not all frames got a plane fitted, as shown by %s',...
            symbolByPlane(2));
        %plot different symbol for planes that don't have a plane fitted
        scatter(t(~planeFitted),nSpots(~planeFitted,i),200)
%            repmat(coloursByStage(i,:),sum(planeFitted),1), symbolByPlane(2));
    else
        error('wrong number of fitted planes compared to frames');
    end
end
legend(stagesLegend,'Location','southwest');
xlabel('Time (s)');
ylabel('Number of spots');
title('Performance of each stage of KiT in tracking kinetochores');
grid on;
box on;
set(gca,'fontsize',20);
savename =  sprintf('%s/%sNumSpotsDiagnosticPlot',job.movieDirectory,...
    job.ROI(1).movie);
savename(regexp(savename,'[.]'))=[];
print(sprintf('%s.eps',savename),'-depsc');

%%%%%%%%%%%%
%call other diagnostics too
if isfield(job,'dataStruct') %these are meaningless if didnt manage this
kitCheckTracks(job);
kitPlotHowManyTracksOfGivenLength(job, channel);
kitMoreDiagnosticPlots(job,channel,1);
kitVisualiseTrackedAndGroupedSpots(job,movie,channel,...
    round(job.metadata.nFrames/2));
end

if applyTests
%% Test1: good number of spots
avSpotsByStage = mean(nSpots,1);
assert(all(avSpotsByStage)>30, 'Want at least 15 sister pairs');

%% Test2: consistency between stages
assert(all(avSpotsByStage(2:nStages)./avSpotsByStage(1:(nStages-1)))<2,...
    'Should not suddenly double number of spots');
assert(all(avSpotsByStage(2:nStages)./avSpotsByStage(1:(nStages-1)))...
    >0.7, 'Should not suddenly lose lots of spots');
%% Test3: consistency across the movie
maxSpotsByStage = max(nSpots,[],1);
assert(all((maxSpotsByStage-avSpotsByStage)./avSpotsByStage<1),...
    ['Check if the movie bleachs badly over time,'...
    'since detections are inconsistent over time']);
end
end
