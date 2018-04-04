function [job,userStatus] = pairTracks(job,opts)
% PAIRTRACKS Manual pairing of tracks in 3D time series.
%
% Copyright (c) 2017 C. A. Smith

% predefine track length to show
dragTails = 5;

%% GET REQUIRED IMAGE AND METADATA
[md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie),job.metadata);
movieIdx = job.index;

% get crop information, if any
crop = job.ROI(movieIdx).crop;
cropSize = job.ROI(movieIdx).cropSize;
if isempty(crop)
    cropSize = md.frameSize;
end
% specify RGB channel order [R G B]
chanOrder = opts.chanOrder;
% calculate size of cropped region in pixels based on maxSisSep
pixelSize = job.metadata.pixelSize(1:3);
cropRange = 1.25*repmat(opts.maxSisSep,1,3)./pixelSize;
cropRange = round(cropRange);
chrShift = job.options.chrShift.result;
% find chrShift rounded to nearest pixel
pixChrShift = cellfun(@times,chrShift,repmat({[pixelSize pixelSize]},3),'UniformOutput',0);
nFrames = job.metadata.nFrames;

%% GET IMAGE AND COORDINATE INFORMATION

% get trackList
% check whether or not this movie has an initCoord
if ~isfield(job.dataStruct{opts.plotChan},'initCoord')
  kitLog('No initCoord present. Skipping movie.');
  trackList = [];
% check whether or not this movie has an initCoord
elseif isfield(job.dataStruct{opts.plotChan},'failed') && job.dataStruct{opts.plotChan}.failed
  kitLog('Tracking failed to produce an initCoord. Skipping movie.');
  trackList = [];
  % check whether or not this movie has any tracks
elseif ~isfield(job.dataStruct{opts.plotChan},'trackList') || isempty(job.dataStruct{opts.plotChan}.trackList(1).coords);
  kitLog('No trackList present. Skipping movie.');
  trackList = [];
else
  initCoord = job.dataStruct{opts.plotChan}.initCoord;  
  trackList = job.dataStruct{opts.plotChan}.trackList;
  tracks = job.dataStruct{opts.plotChan}.tracks;
end
nTracks = length(trackList);

% check that redo hasn't been incorrectly provided here
if ~isfield(job.dataStruct{opts.plotChan},'sisterList') || ...
        isempty(job.dataStruct{opts.plotChan}.sisterList(1).trackPairs)
  opts.redo = 1;
end
% pre-allocate index vectors
if opts.redo
  pairedIdx = [];
  sisterIdxArray = [];
  defaultIdx = 1:nTracks;
else
  % get list of featIdx from current sisterList
  sisterIdxArray = job.dataStruct{opts.plotChan}.sisterList(1).trackPairs(:,1:2);
  pairedIdx = sisterIdxArray(:)';
  defaultIdx = setdiff(1:nTracks,pairedIdx(:));
end
unpairedIdx = [];
doublePairingIdx = [];

%% PLOT IMAGE, REQUEST INPUT

if ~isempty(trackList)
  % produce image file
  fullImg = zeros([cropSize(1:2),3]);
  for iChan = opts.imageChans
    img(:,:,:,iChan) = kitReadImageStack(reader, md, 1, iChan, crop, 0);
    fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
    irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
    fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
  end

  % produce figure environment
  figure(1)
  clf
  plotTitle = sprintf('Check cell orientation\nPress: y to continue, n to skip, q to quit.');
  if length(opts.imageChans) == 1
    imshow(fullImg(:,:,chanOrder(opts.imageChans)));
    title(plotTitle,'FontSize',14)
  else
    imshow(fullImg)
    title(plotTitle,'FontSize',14)
  end
  [~,~,buttonPress] = ginput(1);

  if buttonPress == 110 %110 corresponds to user pressing n key
    unallocatedIdx = [];
  elseif buttonPress == 121 %121 corresponds to user pressing y key
    unallocatedIdx = defaultIdx;
    kitLog('Cell accepted. Continuing with pairing of sisters.');
  elseif buttonPress == 113 %113 corresponds to user pressing q key
    userStatus = 'userPaused';
    return
  else
    unallocatedIdx = defaultIdx;
    kitLog('Another key other than y, n or q pressed. Continuing with kinetochore pairing.');
  end
  % preconfigure progress information
  prog = kitProgress(0);
  
else
  % if not coords, make unallocatedIdx empty to avoid running spot finding
  unallocatedIdx = [];

end

plotTitle = sprintf('Locate the white track''s sister.\nLeft click on an: unallocated (g), ignored (y) or pre-paired (r) track.\nPress backspace in there is no clear sister track.');

while ~isempty(unallocatedIdx)
    
    % give progress information
    nRemaining = size(unallocatedIdx,2);
    prog = kitProgress((nTracks-nRemaining)/nTracks,prog);
    
    % take the earliest trackIdx
    iTrack = unallocatedIdx(1);
    
    % get track information and check whether it is long enough
    thisTrack = tracks(iTrack);
    if diff(thisTrack.seqOfEvents(:,1)) < (2*dragTails+1)
      unallocatedIdx = setdiff(unallocatedIdx,iTrack);
      unpairedIdx = [unpairedIdx iTrack];
      continue
    end
      
    % find timepoints at and around centre of the track being paired
    tP = ceil(mean(thisTrack.seqOfEvents(:,1)));
    bef = max(1,tP-dragTails);
    aft = min(nFrames,tP+dragTails);
    
    % get nnDist information
    unallNNidx    = getNNdistIdx(trackList,iTrack,unallocatedIdx,opts.maxSisSep,bef:aft);
    pairedNNidx   = getNNdistIdx(trackList,iTrack,pairedIdx,opts.maxSisSep,bef:aft);
    unpairedNNidx = getNNdistIdx(trackList,iTrack,unpairedIdx,opts.maxSisSep,bef:aft);
    
    % if spot has no nearest neighbour, remove and continue
    if isempty([unallNNidx pairedNNidx unpairedNNidx])
        kitLog('No spots located within %.1fµm of track %i. Track not paired.',opts.maxSisSep,iTrack)
        unpairedIdx = [unpairedIdx iTrack];
        unallocatedIdx(1) = [];
        continue
    end
    
    % produce image file
    fullImg = zeros([cropSize(1:2),3]);
    for iChan = opts.imageChans
      img(:,:,:,iChan) = kitReadImageStack(reader, md, tP, iChan, crop, 0);
      fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
      irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
      fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
    end
    
    % get origin pixel-coordinates
    startTime = thisTrack.seqOfEvents(1,1);
    featIndx = thisTrack.tracksFeatIndxCG(tP-startTime+1);
    try
    centreCoords = initCoord(tP).allCoordPix(featIndx,1:3);
    catch
        qqq=1;
    end
    centreCoords = round(centreCoords);
    % get pixel-coordinates for each track using nnDistIdx
    originCoordsPix = getPixCoords(tracks,initCoord,iTrack,bef:aft);
    unallCoordsPix = getPixCoords(tracks,initCoord,unallNNidx,bef:aft);
    pairedCoordsPix = getPixCoords(tracks,initCoord,pairedNNidx,bef:aft);
    unpairedCoordsPix = getPixCoords(tracks,initCoord,unpairedNNidx,bef:aft);
    
    % concatenate all coordinates
    
    
    % compile all coords for later - NOT SURE HOW TO HANDLE THIS AT THE MOMENT
%     frameCoordsPix = [originCoordsPix; unallCoordsPix; pairedCoordsPix; unpairedCoordsPix];
    
    % IMAGE PROCESSING
    % predesignate cropImg structure
    coordRange = [max(1,centreCoords(1)-cropRange(1)) min(centreCoords(1)+cropRange(1),cropSize(2));...
                  max(1,centreCoords(2)-cropRange(2)) min(centreCoords(2)+cropRange(2),cropSize(1));...
                  max(1,centreCoords(3)-opts.zProjRange) min(centreCoords(3)+opts.zProjRange,cropSize(3))];
    cropImg = zeros(coordRange(2,2)-coordRange(2,1)+1,...
                    coordRange(1,2)-coordRange(1,1)+1,...
                    3);
    
	% plot image in figure 1
    figure(1)
    clf
                
    switch opts.mode
        
        case 'zoom'
    
            % get zoomed image centred at origin coordinate
            tempImg = img(coordRange(2,1):coordRange(2,2), coordRange(1,1):coordRange(1,2), :, :);
            for iChan = opts.imageChans
                if iChan ~= opts.plotChan
                    iCoordRange = coordRange - pixChrShift{opts.plotChan,iChan}(1:2);
                else
                    iCoordRange = coordRange;
                end
                cropImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, iCoordRange(3,1):iCoordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(cropImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                cropImg(:,:,chanOrder(iChan)) = imadjust(cropImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end

            if length(opts.imageChans) == 1
                imshow(cropImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(cropImg)
            end
            title(plotTitle,'FontSize',14)

            % PLOTTING ZOOMED COORDINATES
            hold on
            % origin coordinates in white
            plot(originCoordsPix{1}(:,1) - coordRange(1,1)+1, originCoordsPix{1}(:,2) - coordRange(2,1)+1,'w','LineWidth',2);
            % unallocated in green
            for iIdx = 1:length(unallCoordsPix)
              plot(unallCoordsPix{iIdx}(:,1) - coordRange(1,1)+1, unallCoordsPix{iIdx}(:,2) - coordRange(2,1)+1,'g','LineWidth',1);
            end
            % unpaired in yellow
            for iIdx = 1:length(unpairedCoordsPix)
              plot(unpairedCoordsPix{iIdx}(:,1) - coordRange(1,1)+1, unpairedCoordsPix{iIdx}(:,2) - coordRange(2,1)+1,'y','LineWidth',1);
            end
            % paired in red
            for iIdx = 1:length(pairedCoordsPix)
              plot(pairedCoordsPix{iIdx}(:,1) - coordRange(1,1)+1, pairedCoordsPix{iIdx}(:,2) - coordRange(2,1)+1,'r','LineWidth',1);
            end
        
        case 'full'
        
            get image
            tempImg = img;
            for iChan = opts.imageChans
                fullImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, coordRange(3,1):coordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end

            if length(opts.imageChans) == 1
                imshow(fullImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(fullImg)
            end
            title(plotTitle,'FontSize',14)
        
            % PLOTTING COORDINATES
            hold on
            % origin coordinates in white
            plot(originCoordsPix{1}(:,1), originCoordsPix{1}(:,2),'w','LineWidth',2);
            % unallocated in green
            for iIdx = 1:length(unallCoordsPix)
              plot(unallCoordsPix{iIdx}(:,1), unallCoordsPix{iIdx}(:,2),'g','LineWidth',1);
            end
            % unpaired in yellow
            for iIdx = 1:length(unpairedCoordsPix)
              plot(unpairedCoordsPix{iIdx}(:,1), unpairedCoordsPix{iIdx}(:,2),'y','LineWidth',1);
            end
            % paired in red
            for iIdx = 1:length(pairedCoordsPix)
              plot(pairedCoordsPix{iIdx}(:,1), pairedCoordsPix{iIdx}(:,2),'r','LineWidth',1);
            end
        
        case 'dual'
            
            subplot(2,4,[3,4,7,8])
            % get full image
            tempImg = img;
            for iChan = opts.imageChans
                fullImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, coordRange(3,1):coordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end
            
            if length(opts.imageChans) == 1
                imshow(fullImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(fullImg)
            end
            
            % PLOTTING COORDINATES
            hold on
            % origin coordinates in white
            plot(originCoordsPix{1}(:,1), originCoordsPix{1}(:,2),'w','LineWidth',2);
            % unallocated in green
            for iIdx = 1:length(unallCoordsPix)
              plot(unallCoordsPix{iIdx}(:,1), unallCoordsPix{iIdx}(:,2),'g','LineWidth',0.5);
            end
            % unpaired in yellow
            for iIdx = 1:length(unpairedCoordsPix)
              plot(unpairedCoordsPix{iIdx}(:,1), unpairedCoordsPix{iIdx}(:,2),'y','LineWidth',0.5);
            end
            % paired in red
            for iIdx = 1:length(pairedCoordsPix)
              plot(pairedCoordsPix{iIdx}(:,1), pairedCoordsPix{iIdx}(:,2),'r','LineWidth',0.5);
            end
            
            subplot(2,4,[1,2,5,6])
            % get zoomed image
            tempImg = img(coordRange(2,1):coordRange(2,2), coordRange(1,1):coordRange(1,2), :, :);
            for iChan = opts.imageChans
                cropImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, coordRange(3,1):coordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(cropImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                cropImg(:,:,chanOrder(iChan)) = imadjust(cropImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end
            
            if length(opts.imageChans) == 1
                imshow(cropImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(cropImg)
            end
            title(plotTitle,'FontSize',14)
            
            % PLOTTING ZOOMED COORDINATES
            hold on
            % origin coordinates in white
            plot(originCoordsPix{1}(:,1) - coordRange(1,1)+1, originCoordsPix{1}(:,2) - coordRange(2,1)+1,'w','LineWidth',2);
            % unallocated in green
            for iIdx = 1:length(unallCoordsPix)
              plot(unallCoordsPix{iIdx}(:,1) - coordRange(1,1)+1, unallCoordsPix{iIdx}(:,2) - coordRange(2,1)+1,'g','LineWidth',0.5);
            end
            % unpaired in yellow
            for iIdx = 1:length(unpairedCoordsPix)
              plot(unpairedCoordsPix{iIdx}(:,1) - coordRange(1,1)+1, unpairedCoordsPix{iIdx}(:,2) - coordRange(2,1)+1,'y','LineWidth',0.5);
            end
            % paired in red
            for iIdx = 1:length(pairedCoordsPix)
              plot(pairedCoordsPix{iIdx}(:,1) - coordRange(1,1)+1, pairedCoordsPix{iIdx}(:,2) - coordRange(2,1)+1,'r','LineWidth',0.5);
            end
            
    end
    
    
    % USER INPUT AND SISTER-SISTER PAIRING
    
    % uses a crosshair to allow user to provide a pair of coordinates
    [userX,userY,key] = ginput(1);
    
    if ismember(key,1:3)
        
        % correct user-provided information if zoomed
        if any(strcmp(opts.mode,{'dual','zoom'}))
            userX = userX + coordRange(1,1)-1; userY = userY + coordRange(2,1)-1;
        end
        
        % find the user-selected track
        % distance relative to origin track
        originUserNNdist = createDistanceMatrix([userX userY],originCoordsPix{1}(:,1:2));
        % unallocated
        tempDist = cat(1,unallCoordsPix{:});
        if ~isempty(tempDist)
          unallUserNNdist = createDistanceMatrix([userX userY],tempDist(:,1:2));
        else
          unallUserNNdist = cropSize(1);
        end
        % paired
        tempDist = cat(1,pairedCoordsPix{:});
        if ~isempty(tempDist)
          pairedUserNNdist = createDistanceMatrix([userX userY],tempDist(:,1:2));
        else
          pairedUserNNdist = cropSize(1); % maximum distance two tracks can be apart
        end
        % and unpaired
        tempDist = cat(1,unpairedCoordsPix{:});
        if ~isempty(tempDist)
          unpairedUserNNdist = createDistanceMatrix([userX userY],tempDist(:,1:2));
        else
          unpairedUserNNdist = cropSize(1);
        end
        
        % find which distribution holds the minimum
        testDist = originUserNNdist(:);
        minDist = 'orig';
        if min(testDist(:)) > min(unallUserNNdist(:));
          testDist = unallUserNNdist;
          minDist = 'unal';
        end
        if min(testDist(:)) > min(unpairedUserNNdist(:));
          testDist = unpairedUserNNdist;
          minDist = 'unpa';
        end
        if min(testDist(:)) > min(pairedUserNNdist(:));
          minDist = 'pair';
        end
        
        % if left click, add; if right click, remove
        if key == 3
            
            kitLog('Sorry, track removal is not currently implemented. Please bear with us.')
            continue
            
        else
            % get the trackID for the nearest track
            switch minDist

                % check user hasn't given the original cross
                case 'orig'
                    kitLog('Either the white track was selected, or the selected location was equidistant from two tracks. Please try again.')
                    continue

                case 'unal'

                    % find the track containing the closest point
                    nnIdx = find(unallUserNNdist == min(unallUserNNdist(:)));
                    nnIdx = floor(nnIdx/(2*dragTails+1)) + 1;
                    nnTrackIdx = unallNNidx(nnIdx);

                    % plot a thick line over the selected track
                    plot(unallCoordsPix{nnIdx}(:,1) - coordRange(1,1)+1, unallCoordsPix{nnIdx}(:,2) - coordRange(2,1)+1,'g','LineWidth',2);

                    % remove new coordinates from unallocated lists, add to paired
                    unallocatedIdx = setdiff(unallocatedIdx,[iTrack nnTrackIdx]);
                    pairedIdx = [pairedIdx iTrack nnTrackIdx];

                case 'unpa'

                    % find the track containing the closest point
                    nnIdx = find(unpairedUserNNdist == min(unpairedUserNNdist(:)));
                    nnIdx = floor(nnIdx/(2*dragTails+1)) + 1;
                    nnTrackIdx = unpairedNNidx(nnIdx);

                    % plot a thick line over the selected track
                    plot(unpairedCoordsPix{nnIdx}(:,1) - coordRange(1,1)+1, unpairedCoordsPix{nnIdx}(:,2) - coordRange(2,1)+1,'y','LineWidth',2);

                    % remove new coordinates from unpaired and unallocated lists,
                    % add to paired
                    unallocatedIdx = setdiff(unallocatedIdx,iTrack);
                    unpairedIdx = setdiff(unpairedIdx,nnTrackIdx);
                    pairedIdx = [pairedIdx iTrack nnTrackIdx];

                case 'pair'

                    % find the track containing the closest point
                    nnIdx = find(pairedUserNNdist == min(pairedUserNNdist(:)));
                    nnIdx = floor(nnIdx/(2*dragTails+1)) + 1;
                    nnTrackIdx = pairedNNidx(nnIdx);

                    % plot a thick line over the selected track
                    plot(pairedCoordsPix{nnIdx}(:,1) - coordRange(1,1)+1, pairedCoordsPix{nnIdx}(:,2) - coordRange(2,1)+1,'r','LineWidth',2);

                    % remove new coordinate from unallocated, add to paired
                    unallocatedIdx = setdiff(unallocatedIdx,iTrack);
                    pairedIdx = [pairedIdx iTrack];

                    % locate the original pair for the new paired index
                    [iRow,iCol] = find(sisterIdxArray == nnTrackIdx);
                    oldPairIdx = sisterIdxArray(iRow,setdiff([1 2],iCol));

                    % add originally-paired index to unallocated, remove old pair
                    % from sister array
                    unallocatedIdx = [unallocatedIdx oldPairIdx];
                    pairedIdx = setdiff(pairedIdx,oldPairIdx);
                    sisterIdxArray(iRow,:) = [];

                    % save that the common sister has been chosen twice - this may
                    % not be necessary, but keep for the time being
                    doublePairingIdx = [doublePairingIdx; oldPairIdx nnTrackIdx];

                otherwise
                    kitLog('Unexpected error occurred. Please try again. If the problem persists, press ctrl+c and request support.')
                    continue

            end

            % save sister pair indices
            kitLog('Tracks paired: %i and %i.',iTrack,nnTrackIdx)
            sisterIdxArray = [sisterIdxArray; iTrack nnTrackIdx];
        end
        
    else
        % case where no spot selected
        unpairedIdx = [unpairedIdx iTrack];
        unallocatedIdx = setdiff(unallocatedIdx,iTrack);
        kitLog('No sister paired with %i.',iTrack)
    end
    
end

prog = kitProgress(1,prog);

userStatus = 'completed';

%% STORE INFORMATION

% set up empty sisterLists
emptySisterList = struct('trackPairs',[],...
                        'coords1',[],...
                        'coords2',[],...
                        'sisterVectors',[],...
                        'distances',[]);
                    
for iChan = opts.coordChans
    
    % construct new sisterList
    sisterList = emptySisterList;
    sisterList(1).trackPairs = sisterIdxArray;
    
    % get microscope coordinates and standard deviations in µm
    for iSis = 1:size(sisterIdxArray,1)
        
        sisterList(iSis).coords1 = trackList(sisterIdxArray(iSis,1)).coords;
        sisterList(iSis).coords2 = trackList(sisterIdxArray(iSis,2)).coords;
        sisterList(iSis).distances = sqrt(sum((sisterList(iSis).coords1 - sisterList(iSis).coords2).^2,2));
        trackList(sisterIdxArray(iSis,1)).sister = sisterIdxArray(iSis,2);
        trackList(sisterIdxArray(iSis,2)).sister = sisterIdxArray(iSis,1);
        
    end
    
    job.dataStruct{iChan}.sisterList = sisterList;
    job.dataStruct{iChan}.trackList = trackList;
    job = kitSaveJob(job);
    
end
end
   
function distIdx = getNNdistIdx(trackList,origIdx,otherIdx,maxSisSep,tPs)

    distIdx = [];
    coords = [];
    for iTime = tPs
        
      for jTrack = 1:length(trackList)
        coords(jTrack,1:3) = trackList(jTrack).coords(iTime,1:3);
      end
    
      nnDists = createDistanceMatrix(coords(origIdx,:),coords(otherIdx,:));
      nnCoords = (nnDists < maxSisSep & nnDists > 0);
      distIdx = [distIdx otherIdx(nnCoords)];
      
    end
    distIdx = unique(distIdx);
    distIdx = distIdx(:)';
    
end

function pixCoords = getPixCoords(tracks,initCoord,allIdx,tPs)
    
    featIndx = nan(1,length(initCoord));
    pixCoords = repmat({nan(length(tPs),3)},length(allIdx),1);
    for idx = 1:length(allIdx)
      
      startTime = tracks(allIdx(idx)).seqOfEvents(1,1);
      endTime = tracks(allIdx(idx)).seqOfEvents(2,1);
      featIndx(startTime:endTime) = tracks(allIdx(idx)).tracksFeatIndxCG;
      
      c = 1;
      for iTime = tPs
        if isnan(featIndx(iTime)) || featIndx(iTime) == 0
          pixCoords{idx}(c,:) = nan;
        else
            try
          pixCoords{idx}(c,:) = initCoord(iTime).allCoordPix(featIndx(iTime),1:3);
            catch
                qqq=1;
            end
        end
        c = c+1;
      end
    end

end
    
    