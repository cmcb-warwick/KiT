function [job,userStatus] = pairSpots(job,opts)
% PAIRSPOTS Manual pairing of spots in single timepoint z-stacks.
%
% Copyright (c) 2018 C. A. Smith

% Check form of ROIs
if length(job.ROI)>1
    job.ROI = job.ROI(job.index);
end

%% GET REQUIRED IMAGE AND METADATA
[md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);

% get channel information
for iCh = 1:length(job.dataStruct)
    opts.coordChans(iCh) = isfield(job.dataStruct{iCh},'initCoord');
end
opts.coordChans = find(opts.coordChans);
% get crop information, if any
crop = job.ROI.crop;
cropSize = job.ROI.cropSize;
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
pixChrShift = cellfun(@times,chrShift,repmat({[pixelSize pixelSize]},size(chrShift,1)),'UniformOutput',0);

%% GET IMAGE AND COORDINATE INFORMATION

% check whether cell has been filtered out by user
if isfield(job,'keep') && ~job.keep
  kitLog('User selected this movie to be ignored in analysis. Skipping movie.');
  coords = [];
  coordsPix = [];
  skip = 1;
% get coordinates in both µm and pixels
% check whether or not this movie has a dataStruct
elseif ~isfield(job,'dataStruct')
  kitLog('No dataStruct present. Skipping movie.');
  coords = [];
  coordsPix = [];
  skip = 1;
% check whether or not this movie has an initCoord
elseif ~isfield(job.dataStruct{opts.plotChan},'initCoord')
  kitLog('No initCoord present. Skipping movie.');
  coords = [];
  coordsPix = [];
  skip = 1;
% check whether or not this movie has an initCoord
elseif isfield(job.dataStruct{opts.plotChan},'failed') && job.dataStruct{opts.plotChan}.failed
  kitLog('Tracking failed to produce an initCoord. Skipping movie.');
  coords = [];
  coordsPix = [];
  skip = 1;
else
    isTransposed = any(size(kitReadImageStack(reader,md,1,2,crop,0)) ~= cropSize);
    if isTransposed
        xyz_ind = [1 2 3];
    else
        xyz_ind = [2 1 3];
    end
    coords = job.dataStruct{opts.plotChan}.initCoord(1).allCoord(:,xyz_ind);
    coordsPix = job.dataStruct{opts.plotChan}.initCoord(1).allCoordPix(:,xyz_ind);
    skip = 0;
end
nCoords = size(coords,1);

%% PLOT IMAGE, REQUEST INPUT

if ~isempty(coords)
    
  % check that redo hasn't been incorrectly provided here
  if ~isfield(job.dataStruct{opts.plotChan},'sisterList') || ...
          isempty(job.dataStruct{opts.plotChan}.sisterList(1).trackPairs)
    opts.redo = 1;
  end
  % pre-allocate index vectors
  if opts.redo
    pairedIdx = [];
    sisterIdxArray = [];
    unallocatedIdx = 1:nCoords;
  else
    % get list of featIdx from current sisterList
    sisterIdxArray = job.dataStruct{opts.plotChan}.sisterList(1).trackPairs(:,1:2);
    pairedIdx = sisterIdxArray(:)';
    unallocatedIdx = setdiff(1:nCoords,pairedIdx(:));
  end
  unpairedIdx = [];
  doublePairingIdx = [];

try  
  % produce image file
  fullImg = zeros([cropSize([1 2]) 3]);
  for iChan = opts.imageChans
    img(:,:,:,iChan) = kitReadImageStack(reader, md, 1, iChan, crop, 0);
    fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
    irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
    fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
  end
  isTransposed = 0; %use this to check later which case we are in
catch ME %image may be transposed
  % produce image file
  fullImg = zeros([cropSize([2 1]) 3]);
  for iChan = opts.imageChans
    img(:,:,:,iChan) = kitReadImageStack(reader, md, 1, iChan, crop, 0);
    fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
    irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
    fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
  end
  isTransposed = 1; %use this to check later
end
  
  % produce figure environment
  f = figure(1);
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
  
  if exist('f','var'); close(f); end

  if buttonPress == 110 %110 corresponds to user pressing n key
    unallocatedIdx = [];
  elseif buttonPress == 121 %121 corresponds to user pressing y key
    kitLog('Cell accepted. Continuing with pairing of sisters.');
  elseif buttonPress == 113 %113 corresponds to user pressing q key
    userStatus = 'userPaused';
    return
  else
    kitLog('Another key other than y, n or q pressed. Continuing with kinetochore pairing.');
  end
  
  % preconfigure progress information
  prog = kitProgress(0);
  
else
  % if not coords, make unallocatedIdx empty to avoid running spot finding
  unallocatedIdx = [];
  sisterIdxArray = [];
end

plotTitle = sprintf('Locate the white cross'' sister.\nClick on an: unallocated (g), ignored (y) or pre-paired (r) cross.\nPress backspace if none applicable.');

while ~isempty(unallocatedIdx)
    
    % give progress information
    nRemaining = length(unallocatedIdx);
    prog = kitProgress((nCoords-nRemaining)/nCoords,prog);
    
    % take the earliest spotIdx
    iCoord = unallocatedIdx(1);

    % get nnDist information
    unallNNidx    = getNNdistIdx(coords,iCoord,unallocatedIdx,opts.maxSisSep);
    pairedNNidx   = getNNdistIdx(coords,iCoord,pairedIdx,opts.maxSisSep);
    unpairedNNidx = getNNdistIdx(coords,iCoord,unpairedIdx,opts.maxSisSep);
    
    % if spot has no nearest neighbour, remove and continue
    if isempty([unallNNidx pairedNNidx unpairedNNidx])
        kitLog('No spots located within %.1fµm of spot %i. Spot not paired.',opts.maxSisSep,iCoord)
        unpairedIdx = [unpairedIdx iCoord];
        unallocatedIdx(1) = [];
        continue
    end
    
    % get origin pixel-coordinates
    iCoordsPix = coordsPix(iCoord,:);
    centreCoords = round(iCoordsPix);
    % get pixel-coordinates for each nnDistIdx
    unallCoordsPix = coordsPix(unallNNidx,:);
    pairedCoordsPix = coordsPix(pairedNNidx,:);
    unpairedCoordsPix = coordsPix(unpairedNNidx,:);
    % compile all coords for later
    frameCoordsPix = [iCoordsPix; unallCoordsPix; pairedCoordsPix; unpairedCoordsPix];
    
    % IMAGE PROCESSING
    % predesignate cropImg structure
    coordRange = [max(1,centreCoords(1)-cropRange(1)) min(centreCoords(1)+cropRange(1),cropSize(1));...
                  max(1,centreCoords(2)-cropRange(2)) min(centreCoords(2)+cropRange(2),cropSize(2));...
                  max(1,centreCoords(3)-opts.zProjRange) min(centreCoords(3)+opts.zProjRange,cropSize(3))];
    if isTransposed
          cropImg = zeros(coordRange(2,2)-coordRange(2,1)+1,...
                coordRange(1,2)-coordRange(1,1)+1,...
                3);
    else
          cropImg = zeros(coordRange(1,2)-coordRange(1,1)+1,...
            coordRange(2,2)-coordRange(2,1)+1,...
            3);
    end
	% plot image in figure 1
    f = figure(1);
    clf
    set(f,'Position',[100 100 800 600]);
                
    switch opts.mode
        
        case 'zoom'
    
            % get zoomed image centred at origin coordinate
            if isTransposed
                tempImg = img(coordRange(2,1):coordRange(2,2), coordRange(1,1):coordRange(1,2), :, :);
            else
                tempImg = img(coordRange(1,1):coordRange(1,2), coordRange(2,1):coordRange(2,2), :, :);                
            end
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
            scatter(iCoordsPix(:,1) - coordRange(2,1)+1, iCoordsPix(:,2) - coordRange(1,1)+1,'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,1) - coordRange(2,1)+1, unallCoordsPix(:,2) - coordRange(1,1)+1,'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,1) - coordRange(2,1)+1, unpairedCoordsPix(:,2) - coordRange(1,1)+1,'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,1) - coordRange(2,1)+1, pairedCoordsPix(:,2) - coordRange(1,1)+1,'xr','sizeData',200,'LineWidth',1.25);
        
        case 'full'
        
            % get image
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
            scatter(iCoordsPix(:,1),iCoordsPix(:,2),'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,1),unallCoordsPix(:,2),'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,1),unpairedCoordsPix(:,2),'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,1),pairedCoordsPix(:,2),'xr','sizeData',200,'LineWidth',1.25);
        
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
            if isTransposed
                scatter(iCoordsPix(:,1),iCoordsPix(:,2),'xw','sizeData',200,'LineWidth',1.25);
            else    
                scatter(iCoordsPix(:,2),iCoordsPix(:,1),'xw','sizeData',200,'LineWidth',1.25);
            end
            
            subplot(2,4,[1,2,5,6])
            % get zoomed image
            if isTransposed
                tempImg = img(coordRange(2,1):coordRange(2,2), coordRange(1,1):coordRange(1,2), :, :);
            else %again due to transpose. Image has been flipped and crops should be too.
                tempImg = img(coordRange(1,1):coordRange(1,2), coordRange(2,1):coordRange(2,2), :, :);
            end
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
            if ~isTransposed
            % origin coordinates in white
            scatter(iCoordsPix(:,2) - iCoordRange(2,1)+1, iCoordsPix(:,1) - iCoordRange(1,1)+1,'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,2) - iCoordRange(2,1)+1, unallCoordsPix(:,1) - iCoordRange(1,1)+1,'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,2) - iCoordRange(2,1)+1, unpairedCoordsPix(:,1) - iCoordRange(1,1)+1,'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,2) - iCoordRange(2,1)+1, pairedCoordsPix(:,1) - iCoordRange(1,1)+1,'xr','sizeData',200,'LineWidth',1);
            else
                % origin coordinates in white
                scatter(iCoordsPix(:,1) - iCoordRange(1,1)+1, iCoordsPix(:,2) - iCoordRange(2,1)+1,'xw','sizeData',200,'LineWidth',1.25);
                % unallocated in green
                scatter(unallCoordsPix(:,1) - iCoordRange(1,1)+1, unallCoordsPix(:,2) - iCoordRange(2,1)+1,'xg','sizeData',200,'LineWidth',1.25);
                % unpaired in yellow
                scatter(unpairedCoordsPix(:,1) - iCoordRange(1,1)+1, unpairedCoordsPix(:,2) - iCoordRange(2,1)+1,'xy','sizeData',200,'LineWidth',1.25);
                % paired in red
                scatter(pairedCoordsPix(:,1) - iCoordRange(1,1)+1, pairedCoordsPix(:,2) - iCoordRange(2,1)+1,'xr','sizeData',200,'LineWidth',1);
            end
                
    end
    
    
    % USER INPUT AND SISTER-SISTER PAIRING
    
    % uses a crosshair to allow user to provide a pair of coordinates
    [userX,userY,key] = ginput(1);
    
    
    if ismember(key,1:3)
        
        % correct user-provided information if zoomed
        if any(strcmp(opts.mode,{'dual','zoom'}))
            if isTransposed
                userX = userX + iCoordRange(1,1)-1; userY = userY + iCoordRange(2,1)-1;
            else
                userX = userX + iCoordRange(2,1)-1; userY = userY + iCoordRange(1,1)-1;
            end
        end
        
        % find the user-selected cross
        if isTransposed
            coordsInd = [1 2];
        else
            coordsInd = [2 1];
        end
        userNNdist = createDistanceMatrix([userX userY],frameCoordsPix(:,coordsInd));
        [~,userIdx] = min(userNNdist);
        [userIdx,~] = find((coordsPix - repmat(frameCoordsPix(userIdx,:),size(coordsPix,1),1)) == 0);
        userIdx = userIdx(1);
        
        % check user hasn't given the original cross
        if userIdx == iCoord
            kitLog('White cross cannot be selected. Please try again.')
            continue
        end
            
        % check from which group the new spot is, and process accordingly
        if any(userIdx == unallocatedIdx)
            
            % remove new coordinates from unallocated lists, add to paired
            unallocatedIdx = setdiff(unallocatedIdx,[iCoord userIdx]);
            pairedIdx = [pairedIdx iCoord userIdx];
            
        elseif any(userIdx == unpairedIdx)
            
            % remove new coordinates from unpaired and unallocated lists,
            % add to paired
            unallocatedIdx = setdiff(unallocatedIdx,iCoord);
            unpairedIdx = setdiff(unpairedIdx,userIdx);
            pairedIdx = [pairedIdx iCoord userIdx];
            
        elseif any(userIdx == pairedIdx)
            
            % remove new coordinate from unallocated, add to paired
            unallocatedIdx = setdiff(unallocatedIdx,iCoord);
            pairedIdx = [pairedIdx iCoord];
            
            % locate the original pair for the new paired index
            [iRow,iCol] = find(sisterIdxArray == userIdx);
            oldPairIdx = sisterIdxArray(iRow,setdiff([1 2],iCol));
            
            % add originally-paired index to unallocated, remove old pair
            % from sister array
            unallocatedIdx = [unallocatedIdx oldPairIdx];
            pairedIdx = setdiff(pairedIdx,oldPairIdx);
            sisterIdxArray(iRow,:) = [];
            
            % save that the common sister has been chosen twice - this may
            % not be necessary, but keep for the time being
            doublePairingIdx = [doublePairingIdx; oldPairIdx userIdx];
            
        end
        % save sister pair indices
        kitLog('Sisters paired: %i and %i.',iCoord,userIdx)
        sisterIdxArray = [sisterIdxArray; iCoord userIdx];
    
    else % coordinate not paired with anything
        unpairedIdx = [unpairedIdx iCoord];
        unallocatedIdx = setdiff(unallocatedIdx,iCoord);
        kitLog('No sister paired with %i.',iCoord)
    end
    
end

% close the figure
if exist('f','var') && isa(f,'matlab.ui.Figure') && ishandle(f)
    close(f);
end

userStatus = 'completed';
if skip; return; end

%% STORE INFORMATION

% set up empty sisterLists
emptySisterList = struct('trackPairs',[],...
                        'coords1',[],...
                        'coords2',[],...
                        'sisterVectors',[],...
                        'distances',[]);
emptyTracks = struct('tracksFeatIndxCG',0,...
               'tracksCoordAmpCG',[],...
               'seqOfEvents', [1 1 1; 1 2 1],...
               'coordAmp4Tracking',[]);
                    
for iChan = opts.coordChans
    
    % construct new sisterList
    sisterList = emptySisterList;
    sisterList(1).trackPairs(:,1) = 1:size(sisterIdxArray,1);
    sisterList(1).trackPairs(:,2) = sisterList(1).trackPairs(:,1)+size(sisterIdxArray,1);
    
    % get microscope coordinates and standard deviations in µm
    if isfield(job.dataStruct{iChan},'initCoord')
      coords = job.dataStruct{iChan}.initCoord(1).allCoord;
      amps = job.dataStruct{iChan}.initCoord(1).amp;
    else
      coords = [];
      amps = [];
    end
    % get plane coordinates if applicable
    if isfield(job.dataStruct{iChan},'planeFit') && ~isempty(job.dataStruct{iChan}.planeFit(1).plane)
      rotCoords = job.dataStruct{iChan}.planeFit(1).rotatedCoord;
    else
      rotCoords = coords;
    end
    
    if ~isempty(rotCoords)
      for iSis = 1:size(sisterIdxArray,1)
        sisterList(iSis).coords1 = rotCoords(sisterIdxArray(iSis,1),:);
        sisterList(iSis).coords2 = rotCoords(sisterIdxArray(iSis,2),:);
        sisterList(iSis).distances(:,1) = sqrt(sum((sisterList(iSis).coords1(:,1:3) - sisterList(iSis).coords2(:,1:3)).^2,2));
        sisterList(iSis).distances(:,2) = sqrt(sum((sisterList(iSis).coords1(:,4:6) - sisterList(iSis).coords2(:,4:6)).^2,2));
      end
    else
      for iSis = 1:size(sisterIdxArray,1)
        sisterList(iSis).coords1 = NaN(1,6);
        sisterList(iSis).coords2 = NaN(1,6);
        sisterList(iSis).distances = NaN(1,2);
      end
    end
    
    if isempty(sisterIdxArray)
      % assign in the case no sisters found
      tracks = emptyTracks;
    end
    
    % construct new tracks
    if ~isempty(rotCoords)
      for iTrack = 1:length(sisterIdxArray(:))
        
        % ensure tracks are re-allocated for new cell
        tracks(iTrack) = emptyTracks;
        tracks(iTrack).tracksFeatIndxCG = sisterIdxArray(iTrack);
        tracks(iTrack).tracksCoordAmpCG(:,1:3) = coords(sisterIdxArray(iTrack),1:3);
        tracks(iTrack).tracksCoordAmpCG(:,[4 8]) = amps(sisterIdxArray(iTrack),1:2);
        tracks(iTrack).tracksCoordAmpCG(:,5:7) = coords(sisterIdxArray(iTrack),4:6);
        tracks(iTrack).coordAmp4Tracking(:,1:3) = rotCoords(sisterIdxArray(iTrack),1:3);
        tracks(iTrack).coordAmp4Tracking(:,[4 8]) = amps(sisterIdxArray(iTrack),1:2);
        tracks(iTrack).coordAmp4Tracking(:,5:7) = rotCoords(sisterIdxArray(iTrack),4:6);
        
      end
    else
      for iTrack = 1:length(sisterIdxArray(:))
        
        % ensure tracks are re-allocated for new cell
        tracks(iTrack) = emptyTracks;
        tracks(iTrack).tracksFeatIndxCG = sisterIdxArray(iTrack);
        tracks(iTrack).tracksCoordAmpCG = NaN(1,8);
        tracks(iTrack).coordAmp4Tracking = NaN(1,8);
        
      end
    end
                 
    job.dataStruct{iChan}.tracks = tracks;
    job.dataStruct{iChan}.sisterList = sisterList;
    job = kitExtractTracks(job,iChan);
    
end

% store information of which channel was used for pairing
job.options.manualPairChannel = opts.plotChan;
job = kitSaveJob(job);

end
   
function distIdx = getNNdistIdx(coords,origIdx,otherIdx,maxSisSep)

    otherIdx = setdiff(otherIdx,origIdx);
    nnDists = createDistanceMatrix(coords(origIdx,:),coords(otherIdx,:));
    nnCoords = (nnDists < maxSisSep);
    distIdx = otherIdx(nnCoords);
    
end
