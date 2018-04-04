function [jobset,jobs] = chrsChangeChromaticShift(jobset,newShift,varargin)
%CHRSCHANGECHROMATICSHIFT Changes the chromatic shift measurements given to
%    a certain jobset.
%
%    [JOBSET,JOBS] = CHRSCHANGECHROMATICSHIFT(JOBSET,NEWSHIFT,...) Gives a
%    new JOBSET, and JOBS if requested by the user, containing the
%    chromatic shift measurements given in NEWSHIFT.
%
%    Options, defaults in {}:-
%
%    chanVect: {[1 2]} or two numbers between 1 and 3. The two channels for
%         which the new chromatic shift measurement applies, directed from
%         the first to the second number of the vector.
%
%    andJobs: {0} or 1. Whether or not to correct all jobs already
%         processed.
%
%
% Copyright (c) 2017 C. A. Smith

% default options
opts.chanVect = [];
opts.andJobs = 0;
% process user-defined options
opts = processOptions(opts,varargin{:});

%% Change jobset information

% check whether the chromatic shift provided is in cell form, otherwise
% ensure a channel vector is provided
if isempty(opts.chanVect)
    if ~iscell(newShift)
        error('If providing only one chromatic shift vector, need to specify channels.')
    else
        chrShift = newShift;
        source = jobset.options.chrShift.jobset;
        for fromChan=1:2
          for toChan=2:3
            if fromChan<toChan && sum(chrShift{fromChan,toChan})~=0
              % want to ask user to select the chromatic shift-derived jobset
              kitLog('Please choose a jobset pertaining to the chromatic shift between channels %i and %i:',fromChan,toChan);
              [filename,~] = uigetfile('*.mat','Select jobset file');
              kitLog('%s selected.',filename);
              source{fromChan,toChan} = filename;
              source{toChan,fromChan} = filename;
            else
              source{fromChan,toChan} = [];
              source{toChan,fromChan} = [];
            end
          end
        end
    end
    
% check whether the channel vector includes only allowed channels
elseif any(~ismember(opts.chanVect,1:3))
    error('Can only choose channels from 1 to 3.')
    
else % if channel vector provided, proceed as appropriate

    fromChan = opts.chanVect(1);
    toChan = opts.chanVect(2);
    
    % get original chromatic shift cell array
    chrShift = jobset.options.chrShift.result;
    source = jobset.options.chrShift.jobset;
    
    % replace the newShift with only the vector specified in chanVect
    if iscell(newShift)
        warning('Replacing only specified channel from new shift structure.')
        newShift = newShift{fromChan,toChan};
    end
    
    % ensure that only the coordinate measurements are replaced
    lenShift = length(newShift);
    chrShift{fromChan,toChan}(1:lenShift) = newShift;
    
    % replace the mirrored measurement (ie. from toChan to fromChan)
    revShift = newShift;
    revShift(1:3) = revShift(1:3)*-1;
    chrShift{toChan,fromChan}(1:lenShift) = revShift;
    
    % want to ask user to select the chromatic shift-derived jobset
    kitLog('Please choose a jobset pertaining to the chromatic shift between channels %i and %i:',fromChan,toChan);
    [filename,~] = uigetfile('*.mat','Select jobset file');
    kitLog('%s selected.',filename);
    source{fromChan,toChan} = filename;
    source{toChan,fromChan} = filename;
    
end

% change chromatic shift information in the jobset and save
jobset.options.chrShift.result = chrShift;
jobset.options.chrShift.jobset = source;
kitSaveJobset(jobset);

%% Change jobs information if required

if opts.andJobs
    
    % get job structure
    jobs = kitLoadAllJobs(jobset);
    
    % for each job in the jobset, change the chromatic shift information
    for iMov = 1:length(jobs)
        
        % get some basic metadata
        oldChrShift = jobs{iMov}.options.chrShift.result;
        % change the chromatic shift information
        jobs{iMov}.options.chrShift.result = chrShift;
            
        if isempty(opts.chanVect)
            
            for fromChan = 1:2
                for toChan = 2:3
                    
                    % if no change, skip
                    if isequal(chrShift{fromChan,toChan}(:,1:3),oldChrShift{fromChan,toChan}(:,1:3));
                        break
                    end
                    
                    % calculate the difference between new and old chromatic shift measurements
                    diffChrShift = chrShift{fromChan,toChan}(:,1:3) - oldChrShift{fromChan,toChan}(:,1:3);

                    % check for initCoord measures, and change accordingly
                    if ~(length(jobs{iMov}.dataStruct) < toChan) && isfield(jobs{iMov}.dataStruct{toChan},'failed') ...
                            && ~jobs{iMov}.dataStruct{toChan}.failed ...
                            && ~isempty(jobs{iMov}.dataStruct{toChan}.initCoord(1).allCoord)

                        % get dataStruct, and shift the data
                        dataStruct = jobs{iMov}.dataStruct{toChan};
                        dataStruct = shiftData(dataStruct,diffChrShift,jobs{iMov}.metadata);

                    end

                    % save dataStruct
                    jobs{iMov}.dataStruct{toChan} = dataStruct;
                    if isfield(dataStruct,'tracks')
                        jobs{iMov} = kitExtractTracks(jobs{iMov},toChan);
                    end

                end
            end
            
        else
            
            % if no change, skip
            if isequal(chrShift{fromChan,toChan}(:,1:3),oldChrShift{fromChan,toChan}(:,1:3))
                break
            end

            % calculate the difference between new and old chromatic shift measurements
            diffChrShift = chrShift{fromChan,toChan}(:,1:3) - oldChrShift{fromChan,toChan}(:,1:3);
                        
            % check for initCoord measures, and change accordingly
            if ~(length(jobs{iMov}.dataStruct) < toChan) && isfield(jobs{iMov}.dataStruct{toChan},'failed') ...
                    && ~jobs{iMov}.dataStruct{toChan}.failed ...
                    && ~isempty(jobs{iMov}.dataStruct{toChan}.initCoord(1).allCoord)
                
                % get dataStruct, and shift the data
                dataStruct = jobs{iMov}.dataStruct{toChan};
                dataStruct = shiftData(dataStruct,diffChrShift,jobs{iMov}.metadata);
                
            end
            
            % save dataStruct
            jobs{iMov}.dataStruct{toChan} = dataStruct;
            if isfield(dataStruct,'tracks')
                jobs{iMov} = kitExtractTracks(jobs{iMov},toChan);
            end
            
        end
        
        kitSaveJob(jobs{iMov});
        
    end
else
    jobs = NaN;
end

end

%% SUBFUNCTIONS

function dataStruct = shiftData(dataStruct,diffChrShift,metadata)

    % get all dataStruct information
    initCoord = dataStruct.initCoord; % initCoord
    if isfield(dataStruct,'planeFit') % planeFit
        planeFit = dataStruct.planeFit;
    else
        planeFit = [];
    end
    if isfield(dataStruct,'tracks') % tracks
        tracks = dataStruct.tracks;
        nTracks = length(tracks);
    else
        tracks = [];
        nTracks = 0;
    end
    if isfield(dataStruct,'sisterList') && ~isempty(dataStruct.sisterList(1).trackPairs) % sisterList
        sisterList = dataStruct.sisterList;
        nSisters = length(sisterList);
    else
        sisterList = struct('trackPairs',[]);
        nSisters = 0;
    end
    
    % get some metadata
    pixelSize = metadata.pixelSize(1:3);
	nFrames = metadata.nFrames;

    % find which frames have a plane
    framesWiPlane = [];
    for iFrame = 1:nFrames
        if ~isempty(planeFit) && planeFit(iFrame).planeVectorClassifier
            framesWiPlane = [framesWiPlane iFrame];
        end
    end
    % create a rotated chromatic shift structure for later
    rotChrShift = repmat(diffChrShift,nFrames,1);
    % get all coordinate information
    if length(framesWiPlane) > 1

        %get first and last frames to rotate
        firstFrameRotate = framesWiPlane(1);
        lastFrameRotate = framesWiPlane(end);

        %find frames without plane that are before the first frame with plane
        framesBefore1 = 1:firstFrameRotate-1;
        framesAfter1 = lastFrameRotate+1:nFrames;

        %get the coordinate system of each frame with a plane
        coordSystem = zeros(3,3,nFrames);
        coordSystem(:,:,framesWiPlane) = cat(3,planeFit.planeVectors);
        % give any frames without a coordinate system the same as the nearest with a plane
        for iTime = framesBefore1(end:-1:1)
            coordSystem(:,:,iTime) = coordSystem(:,:,iTime+1);
        end
        for iTime = framesAfter1
            coordSystem(:,:,iTime) = coordSystem(:,:,iTime-1);
        end

    else
        % give a standard rotation-free coordinate system
        coordSystem = zeros(3,3,nFrames);
    end

    % shift and rotate initCoord
    for iFrame = 1:nFrames

        % get number of spots, then add the chromatic shift difference to initCoords
        nSpots = initCoord(iFrame).nSpots;
        if nSpots == 0
            continue
        end
        initCoord(iFrame).allCoord(:,1:3) = initCoord(iFrame).allCoord(:,1:3) ...
            + repmat(diffChrShift,nSpots,1);
        initCoord(iFrame).allCoordPix(:,1:3) = initCoord(iFrame).allCoordPix(:,1:3) ...
            + repmat(diffChrShift./pixelSize,nSpots,1);

        % change plane fit
        if sum(coordSystem(:,:,iFrame)==0)==3
            planeFit(iFrame).rotatedCoord = initCoord(iFrame).allCoord;
            planeFit(iFrame).rotatedCoord(:,1:3) = initCoord(iFrame).allCoord(:,1:3) - planeFit(iFrame).planeOrigin(:,1:3);
        else
            rotationMat = inv(coordSystem(:,:,iFrame)); %rotation matrix
            errorPropMat = rotationMat.^2; %error propagation matrix
            planeFit(iFrame).rotatedCoord(:,1:3) = (rotationMat*(initCoord(iFrame).allCoord(:,1:3))')';
            planeFit(iFrame).rotatedCoord(:,4:6) = sqrt((errorPropMat*((initCoord(iFrame).allCoord(:,4:6)).^2)')');
            % form a rotational chromatic shift correction
            rotChrShift(iFrame,1:3) = (rotationMat*(diffChrShift(:,1:3))')';
        end

    end

    % chromatic shift correct tracks
    for iTrack = 1:nTracks
        for iCoord = 1:3
            tracks(iTrack).tracksCoordAmpCG(:,iCoord:8:end) = tracks(iTrack).tracksCoordAmpCG(:,iCoord:8:end) ...
                + diffChrShift(iCoord);
            startTime = tracks(iTrack).seqOfEvents(1,1);
            endTime = tracks(iTrack).seqOfEvents(2,1);
            tracks(iTrack).coordAmp4Tracking(:,iCoord:8:end) = tracks(iTrack).coordAmp4Tracking(:,iCoord:8:end) ...
                + rotChrShift(startTime:endTime,iCoord)';
        end
    end

    % chromatic shift correct sisterList
    for iSis = 1:nSisters
        sisterList(iSis).coords1(:,1:3) = sisterList(iSis).coords1(:,1:3) ...
            + rotChrShift(:,1:3);
        sisterList(iSis).coords2(:,1:3) = sisterList(iSis).coords2(:,1:3) ...
            + rotChrShift(:,1:3);
    end

    % save initCoord, planeFit and tracks information
    dataStruct.initCoord = initCoord;
    dataStruct.planeFit = planeFit;
    if ~isempty(tracks)
        dataStruct.tracks = tracks;
        dataStruct.sisterList = sisterList;
    end

end




