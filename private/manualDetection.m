function spots = manualDetection(movie,metadata,options)
% MANUALDETECTION Allows manual spot-finding in 3D in a given CHANNEL for a
% JOB.
%
% Copyright (c) 2017 C. A. Smith

%% Input + initialization

% get pixel information
pixelSize = metadata.pixelSize(1:3);
is3D = metadata.is3D;
warnDist = options.manualDetect.warningDist;

% movie length
nFrames = metadata.nFrames;
% get image size information
[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);

% get timepoints over which to request user to manually select spots
spacing = options.manualDetect.numFrames-1;
spacing = floor(nFrames/spacing);
timePoints = 1:spacing:nFrames;
if ~ismember(nFrames,timePoints)
  timePoints(end+1) = nFrames;
end

% pre-designate spots structure
spots = cell(nFrames,1);

% turn off warnings off
warning('off','all');
% start progress bar
prog = kitProgress(0);

% loop over the timepoints requested
for iFrame = timePoints

%% GUI to allow user to click on multiple spots: XY

    % get z-projection (if necessary), and contrast
    if is3D
      xyImg = max(movie(:,:,:,iFrame),[],3);
    else
      xyImg = movie(:,:,iFrame);
    end
    irange = stretchlim(xyImg,[0.1 1]);
    xyImg = imadjust(xyImg,irange,[]);
    
    % produce figure environment
    f = figure(1);
    clf
    % plot image
    plotTitle = sprintf('Left-click on centre of all spots\nRight-click on any crosses to delete\nPress d to finish.');
    imshow(xyImg);
    title(plotTitle,'FontSize',14)
    hold on % ensure any coordinates overlaid are indeed overlaid
    set(f,'Position',[100 100 800 600]);
    
    % while loop to allow multiple spots to be added
    ask = 1;
    xySpots = [];
    while ask
      
      % ask user for input
      [userY,userX,key] = ginput(1);
      if isempty(key)
        continue
      end
      
      % check user's input
      switch key
        
        case 100 
          % corresponding to a 'd', break the while loop
          break
          
        case 1
          % get 2D distance to all other spots
          if ~isempty(xySpots)
            dists = abs(xySpots-repmat([userX userY],size(xySpots,1),1));
            dists = dists.*repmat(pixelSize(1:2),size(dists,1),1);
            dists = sqrt(sum(dists.^2,2));
          else
            dists = [];
          end
          
          % check if any distances are within 200nm warning
          if any(dists<warnDist)
            kitLog(['The selected coordinate is within ' num2str(warnDist*1000) 'nm of at least one other. Continue to add this coordinate? (y or n): ']);
            result = input('','s');
          else
            result = 'y';
          end
          
          if strcmp(result,'y')
            % add spot to list
            xySpots = [xySpots; userX userY];
            % draw a green cross where spot is added
            scatter(xySpots(end,2),xySpots(end,1),'gx');
          end
          
        case {2,3}
          
          % check whether there are any spots to delete
          if isempty(xySpots)
            warning('There are no spots to delete. Use the left-click to add spots.')
            continue
          end
          
          % ask user if they are sure about deleting
          kitLog('Are you sure you want to delete this coordinate? (y or n): ');
          result = input('','s');
          if strcmp(result,'y')
            % find closest spot to the coordinates
            dists = abs(xySpots(:,1)-userX);
            [~,r] = min(dists);

            % draw over the deleted coordinate with red cross
            scatter(xySpots(r,2),xySpots(r,1),'rx');

            % remove this coordinate
            xySpots(r,:) = [];
          end
          
        otherwise
              
            continue
          
      end
    
    end
    
    % round coordinates to full pixels
    xySpots = round(xySpots);
    
%% GUI to allow user to click on multiple spots: orthogonal coordinate
    
  if is3D
    
    % get number of spots
    nSpots = size(xySpots,1);
    orthSpots = zeros(nSpots,1);
    % predefine toDelete for non-unique intensities
    toDelete = [];
    for iSpot = 1:nSpots
      
      % calculate pixels needed to get ~500x500nm square around each spot
      r = round(0.25/pixelSize(1));
      xCrop = max(1,xySpots(iSpot,1)-r):min(imageSizeX,xySpots(iSpot,1)+r);
      yCrop = max(1,xySpots(iSpot,2)-r):min(imageSizeY,xySpots(iSpot,2)+r);
      
      % get column of intensities from movie, and find largest sum over xy
      imgCol = movie(xCrop,yCrop,:,iFrame);
      imgCol = sum(imgCol,2); imgCol = sum(imgCol,1);
      maxZslice = find(imgCol == max(imgCol));
      % check that the maximum summed intensity is unique
      if isempty(maxZslice) || length(maxZslice)>1
        toDelete = iSpot;
      else
        orthSpots(iSpot) = maxZslice;
      end
    
    end
    
    % remove problem spots
    xySpots(toDelete,:) = [];
    orthSpots(toDelete,:) = [];
    
    % add good spots to the final list
    spots{iFrame} = [spots{iFrame}; xySpots orthSpots];
    
  else
    
    % if no z-directional information, give only the xy-coordinates
    spots{iFrame} = round(xySpots);
    
  end %is3D
  
  % display progress
  prog = kitProgress(iFrame/nFrames,prog);
  
end

%% Gap filling

% if all timepoints have been processed
if isempty(setdiff(1:nFrames,timePoints))
    return
end

% get gap-filling method
gapMethod = options.manualDetect.gapMethod;

switch gapMethod
    
  case 'linear'
    
    % loop over timepoints (except first)
    for iFrame = [timePoints(1:end-1); timePoints(2:end)]
      
      % get difference between the two
      frameDiff = diff(iFrame);
      
      % get the spots for each frame
      tempSpots1 = spots{iFrame(1)};
      tempSpots2 = spots{iFrame(2)};
      
      % find number of spots possible to add to the list
      if any([isempty(tempSpots1),isempty(tempSpots2)])
        toRun = 0;
      else
        toRun = min(size(tempSpots1,1),size(tempSpots2,1));
      end
      
      while toRun > 0
      
        % find distance between orthogonal coords in adjacent timePoints
        dists = createDistanceMatrix(tempSpots1,tempSpots2);
        dists = abs(dists);
        % find smallest difference between XY and orthogonal spots
        [r,c] = find(dists == min(dists(:)));
        % get coordinate specific differences
        dists = tempSpots2(c,:)-tempSpots1(r,:);
        % get estimated frame-by-frame distance
        frameDisp = dists/frameDiff;
        
        for jFrame = iFrame(1)+1:iFrame(2)-1
          % add these spots to the final list
          spots{jFrame} = [spots{jFrame}; tempSpots1(r,:)+(jFrame-iFrame(1))*frameDisp];
        end
        % delete the coordinates once finished copying
        tempSpots1(r,:) = [];
        tempSpots2(c,:) = [];
        toRun = toRun-1;
      
      end
      
    end

  case 'framewise'
    
    % produce search mask
    r = round(options.maxSearchRadius(1)/(3*pixelSize(1)));
    se = strel('disk', r, 0);
    mask = double(se.getnhood());
    mask(mask == 0) = NaN;
    rz = round(options.maxSearchRadius(1)/(3*pixelSize(3)));
    
    % loop over manual timepoints
    for iFrame = timePoints
      
      % check whether this timepoint is the first or last
      first=(iFrame==1);
      last =(iFrame==nFrames);
      
      % if this frame contains no spots, continue
      if isempty(spots{iFrame})
        continue
      end
      
      % searching ahead
      if ~last
        
        % get list of frames to search
        toCheck = iFrame+1:min(iFrame+ceil(spacing/2),nFrames);
        for jFrame = toCheck
            
          % get source coordinates and num of spots
          startCoords = spots{jFrame-1};
          nSpots = size(startCoords,1);
          for iSpot = 1:nSpots
              
            coords = startCoords(iSpot,:);
            zRange = max(coords(3)-rz,1):min(coords(3)+rz,imageSizeZ);
            
            % make sure that the mask doesn't extend over image boundaries
            if coords(1)-r>0 && coords(1)+r<=imageSizeX && ...
                  coords(2)-r>0 && coords(2)+r<=imageSizeY;
              
              % get original image, and limit to mask
              spotImg = movie(coords(1)-r:coords(1)+r, ...
                    coords(2)-r:coords(2)+r,zRange,jFrame);
              spotImg = max(spotImg,[],3);
              spotImg = mask .* spotImg;
            else
              continue
            end
          
            % find local maximum intensity
            locMaxInt = nanmax(spotImg(:));
            locMax1DIndx = find(spotImg==locMaxInt);
            if isempty(locMax1DIndx)
              continue
            end
            
            for iMax = 1:length(locMax1DIndx)
                
              [locMaxCrd(1),locMaxCrd(2),locMaxCrd(3)] = ind2sub([2*r+1 2*r+1 1],locMax1DIndx(iMax));
              % correct coordinates to full image
              locMaxCrd(1) = locMaxCrd(1)+coords(1)-(r+1);
              locMaxCrd(2) = locMaxCrd(2)+coords(2)-(r+1);
              % find the z-coordinate
              tempCoord = find(movie(locMaxCrd(1),locMaxCrd(2),:,jFrame)==locMaxInt);
              if length(tempCoord)>1 || isnan(tempCoord)
                locMaxCrd(3) = coords(3);
              else
                locMaxCrd(3) = find(movie(locMaxCrd(1),locMaxCrd(2),:,jFrame)==locMaxInt);
              end

              if locMaxCrd(1)>imageSizeX || locMaxCrd(2)>imageSizeY || locMaxCrd(3)>imageSizeZ;
                continue
              end

              % compile coordinates
              spots{jFrame} = [spots{jFrame}; locMaxCrd];
              
            end
          
          end %iSpot
        end %jFrame
      end
      
      % searching before
      if ~first
        
        % get list of frames to search
        toCheck = iFrame-1:-1:iFrame-ceil(spacing/2);
        for jFrame = toCheck
          
          % skip if this timepoint has been previously checked
          if ~isempty(spots{jFrame})
            continue
          end
            
          % get source coordinates and num of spots
          startCoords = spots{jFrame+1};
          nSpots = size(startCoords,1);
          for iSpot = 1:nSpots
            
            coords = startCoords(iSpot,:);
            zRange = max(coords(3)-rz,1):min(coords(3)+rz,imageSizeZ);
          
            % make sure that the mask doesn't extend over image boundaries
            if coords(1)-r>0 && coords(1)+r<=imageSizeX && ...
                  coords(2)-r>0 && coords(2)+r<=imageSizeY;
              % get intensity values of mask-derived spot vicinity
              spotImg = movie(coords(1)-r:coords(1)+r, ...
                    coords(2)-r:coords(2)+r,zRange,jFrame);
              spotImg = max(spotImg,[],3);
              spotImg = mask .* spotImg;
            else
              continue
            end
          
            % find local maximum intensity
            locMaxInt = nanmax(spotImg(:));
            locMax1DIndx = find(spotImg==locMaxInt);
            if isempty(locMax1DIndx)
              continue
            end
              
            for iMax = 1:length(locMax1DIndx)

              [locMaxCrd(1),locMaxCrd(2),locMaxCrd(3)] = ind2sub([2*r+1 2*r+1 1],locMax1DIndx(iMax));
              % correct coordinates to full image
              locMaxCrd(1) = locMaxCrd(1)+coords(1)-(r+1);
              locMaxCrd(2) = locMaxCrd(2)+coords(2)-(r+1);
              tempCoord = find(movie(locMaxCrd(1),locMaxCrd(2),:,jFrame)==locMaxInt);
              if length(tempCoord)>1 || isnan(tempCoord)
                locMaxCrd(3) = coords(3);
              else
                locMaxCrd(3) = find(movie(locMaxCrd(1),locMaxCrd(2),:,jFrame)==locMaxInt);
              end

              if locMaxCrd(1)>imageSizeX || locMaxCrd(2)>imageSizeY || locMaxCrd(3)>imageSizeZ;
                continue
              end

              % compile coordinates
              spots{jFrame} = [spots{jFrame}; locMaxCrd];
              
            end
          
          end %iSpot
        end %jFrame
      end
        
    end
    
  otherwise 

    error('Method of gap-filling not recognised: %s', gapMethod)
      
end

end