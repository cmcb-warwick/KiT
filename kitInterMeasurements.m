function compiledInter = kitInterMeasurements(movies,varargin)
% KITINTERMEASUREMENTS Produces a structure of population-level
% inter-kinetochore measurements over multiple experiments.
%
%    KITINTERMEASUREMENTS(MOVIES,...) Provides all coordinates and inter-
%    kinetochore measurements for all sisters within all cells
%    across all experiments. The resulting structure provides the tools for
%    deriving population-level analyses for an experiment. Options are
%    available.
%
%    Options, defaults in {}:-
%
%    category: {[]} or string. The category from which to calculate intra-
%       measurements.
%
%    channel: 1, 2 or 3. Default is coordinate system channel. The channel 
%       used to make inter-kinetochore measurements.
%
%    prevMeas: {[]} or a structure previously generated. A structure of
%       results from previous experiments to allow new experiment data to
%       be appended.
%
%    spotSelection: {[]} or output from kitSelectData. A structure
%       containing a selection of either sister pair or track IDs per
%       movie, for each experiment. Allows for only specific sisters or
%       spots to be included in the data collection.
%
%
% Copyright (c) 2018 C. A. Smith

% default options
opts.category = [];
opts.channel = [];
opts.prevMeas = [];
opts.spotSelection = [];
% user options
opts = processOptions(opts,varargin{:});

%% Pre-processing input structure

%check structure of movies
if ~iscell(movies{1})
    movies = {movies};
    kitLog('Movie structure provided implies only one experiment. Assuming only one experiment.');
end
%find number of movies
numExpts1 = length(movies);

%check channel information
if isempty(opts.channel)
    ch = movies{1}{1}.options.coordSystemChannel;
else
    ch = opts.channel;
end
    
%process input so that all structs are in cell format
if isempty(opts.spotSelection)
  subset = repmat({[]},numExpts1,1);
  selType = 0;
elseif isstruct(opts.spotSelection) && isfield(opts.spotSelection,'dataType')
  subset = opts.spotSelection.selection;
  switch opts.spotSelection.dataType
    case 'spots' %tracks
      selType = 1;
    case 'sisters'
      selType = 2;
    case 'initCoord'
      selType = 3;
  end
else
  kitLog('Provided spotSelection structure was not derived from kitSelectData. No selection will be imposed.');
  subset = repmat({[]},numExpts1,1);
  selType = 0;
end

%find number of movies and sisters, and ensure they match
numExpts2 = length(subset);
if numExpts1 ~= numExpts2
  error('Have %i spot selections for %i experiments. Please provide spot selection for each experiment.',numExpts2,numExpts1)
end
numExpts = numExpts1;

%% Pre-processing output structure

% assume no pairing until we find one movie that is
paired = 0;
for iExpt = 1:numExpts
    for iMov = 1:length(movies{iExpt})
        paired = isfield(movies{iExpt}{iMov}.dataStruct{ch},'sisterList');
        if paired; break; end
    end
    if paired; break; end
end
if ~paired
    error('Data provided does not contain paired information. Before re-running, run tracking with sister pairing for movies, or use kitManualPairSisters for single-timepoint images.');
end
if selType==3
    subset = opts.spotSelection.selection;
    selType = 1;
end

if isempty(opts.prevMeas)
    
    % make new intra-measurements structure if no previous measurements
    % provided
    allData = kitMakeInterStructure;
    allData = struct2strForm(allData);
    
else
    % get all old data
    allData = struct2strForm(opts.prevMeas);
    
end

% predesignate error arrays
noFail = [];
noSis = [];
noSkip = [];

%% Compiling measurements

for iExpt = 1:numExpts
    
    % get movie and sister list
    theseMovies = movies{iExpt};
    iSubset = subset{iExpt};
    
    for iMov = 1:length(theseMovies)
        
      % get the movie index
      movNum = theseMovies{iMov}.index;
      % check whether there is data in this movie
      if ~isfield(theseMovies{iMov},'dataStruct')
        noFail = [noFail; iExpt movNum];
        continue
      end
      
      % get dataStructs
      dS = theseMovies{iMov}.dataStruct{ch};
      
      % check whether the movie failed
      if (isfield(dS,'failed') && dS.failed) || ~isfield(dS,'initCoord')
        noFail = [noFail; iExpt movNum];
        continue
      end
      
      % check whether the user skipped this movie
      if (isfield(theseMovies{iMov},'keep') && ~theseMovies{iMov}.keep)
        noSkip = [noSkip; iExpt movNum];
        continue
      end
      
      % get basic metadata
      nFrames = theseMovies{iMov}.metadata.nFrames;
      % check if there is a plate fit
      plane = (isfield(dS,'planeFit') && ~isempty(dS.planeFit) && ~isempty(dS.planeFit.planeVectors));
      
      % get initCoord structures - if one doesn't exist trade for an empty
      % initCoord
      iC = dS.initCoord;  
      
      % get spotInts
      if isfield(dS,'spotInt')
        sI = dS.spotInt;
        bg = dS.cellInt.back;
      else
        bg = nan(nFrames,1);
      end
      
      % check whether a sisterList is present, and if it contains any sisters
      if ~isfield(dS,'sisterList') || isempty(dS.sisterList(1).trackPairs)
        noSis = [noSis; iExpt movNum];
        continue
      end
        
      % if no sisters given, go through all sisters in movie
      switch selType
          case 0 % no spot selection
              iSubset = 1:length(dS.sisterList);
          case 1 % using spots/tracks
              iSubset = 1:length(dS.sisterList);
              theseTracks = subset{iExpt}(subset{iExpt}(:,1)==movNum,2)';
          case 2 % using sisters
              iSubset = subset{iExpt}(subset{iExpt}(:,1)==movNum,2)';
      end
        
        for iSis = iSubset
            
          % construct sister pair label
          label = sprintf('%02d%02d%03d',iExpt,iMov,iSis);
            
          % start counter for storing data
          c=1;

          % get sisterLists
          sL = dS.sisterList(iSis);
            
          % get trackID and spotIDs, make spotIDs nan if deselected
          trackIDs = dS.sisterList(1).trackPairs(iSis,1:2);
          spotIDs = nan(nFrames,2);
          for iTrack = 1:2
              if selType~=1 %none or sisters
                  spotIDs(:,iTrack) = dS.trackList(trackIDs(iTrack)).featIndx;
              else %tracks
                  if ismember(trackIDs(iTrack),theseTracks)
                      spotIDs(:,iTrack) = dS.trackList(trackIDs(iTrack)).featIndx;
                  else
                      trackIDs(iTrack) = NaN;
                      switch iTrack
                          case 1
                              sL.coords1(:,1:3) = NaN;
                          case 2
                              sL.coords2(:,1:3) = NaN;
                      end
                  end
              end 
          end
          % filter based on chosen category
          if ~isempty(opts.category) && nFrames==1
              if isfield(theseMovies{iMov},'categories') && ...
                      isfield(theseMovies{iMov}.categories,opts.category)
                  spotIDs = intersect(spotIDs,theseMovies{iMov}.categories.(opts.category));
              else
                  trackIDs = [NaN NaN];
              end
          end
          % if both spots skipped
          if all(isnan(trackIDs))
              continue;
          end

          %% Kinetochore positions
            
          % get microscope coordinates, and plane coordinates if present
          mCoords = nan(nFrames,6);
          pCoords = nan(nFrames,6);
          for iFrame = 1:nFrames
              for iTrk = 1:2
                  % track one stored in (:,1:3), two in (:,4:6)
                  rng = (3*(iTrk-1)+1):3*iTrk;
                  if ~isnan(spotIDs(iFrame,iTrk))
                      mCoords(iFrame,rng) = iC(iFrame).allCoord(spotIDs(iFrame,iTrk),1:3);
                      % check whether or not this movie has a planeFit
                      if plane
                          % rotate coordinates into plane
                          pF = dS.planeFit;
                          if ~isempty(pF(iFrame).planeVectors)
                              coordSystem = pF(iFrame).planeVectors;
                              pCoords(iFrame,:) = mCoords(iFrame,:) - repmat(pF.planeOrigin,1,2);

                              pCoords(iFrame,1:3) = (coordSystem\(pCoords(iFrame,1:3)'))';
                              pCoords(iFrame,4:6) = (coordSystem\(pCoords(iFrame,4:6)'))';
                          end
                      end
                  end
              end
          end
            
          % put data into string format
          newData(c,:) = {'label',label}; c=c+1;
            
          % get microscope coordinates of each spot
          mCoords_x = mCoords(:,[1 4]);
          mCoords_y = mCoords(:,[2 5]);
          mCoords_z = mCoords(:,[3 6]);
          % put data into string format
          newData(c,:) = {'microscope.coords.x',mCoords_x}; c=c+1;
          newData(c,:) = {'microscope.coords.y',mCoords_y}; c=c+1;
          newData(c,:) = {'microscope.coords.z',mCoords_z}; c=c+1;
            
            % get plate coordinates of each spot
            pCoords_x = pCoords(:,[1 4]);
            pCoords_y = pCoords(:,[2 5]);
            pCoords_z = pCoords(:,[3 6]);
            % put data into string format
            newData(c,:) = {'plate.coords.x',pCoords_x}; c=c+1;
            newData(c,:) = {'plate.coords.y',pCoords_y}; c=c+1;
            newData(c,:) = {'plate.coords.z',pCoords_z}; c=c+1;
            
            %% Inter- and intra-kinetochore measurements
            
            micrData = makeMeasurements(mCoords,0);
            % put data into string format
            newData(c,:) = {'microscope.sisSep.x',micrData.sisSep_x}; c=c+1;
            newData(c,:) = {'microscope.sisSep.y',micrData.sisSep_y}; c=c+1;
            newData(c,:) = {'microscope.sisSep.z',micrData.sisSep_z}; c=c+1;
            newData(c,:) = {'microscope.sisSep.twoD',micrData.sisSep_2D}; c=c+1;
            newData(c,:) = {'microscope.sisSep.threeD',micrData.sisSep_3D}; c=c+1;
            
            plateData = makeMeasurements(pCoords,1);
            % put data into string format
            newData(c,:) = {'plate.sisSep.x',plateData.sisSep_x}; c=c+1;
            newData(c,:) = {'plate.sisSep.y',plateData.sisSep_y}; c=c+1;
            newData(c,:) = {'plate.sisSep.z',plateData.sisSep_z}; c=c+1;
            newData(c,:) = {'plate.sisSep.twoD',plateData.sisSep_2D}; c=c+1;
            newData(c,:) = {'plate.sisSep.threeD',plateData.sisSep_3D}; c=c+1;
            newData(c,:) = {'plate.twist.y',plateData.twist_y}; c=c+1;
            newData(c,:) = {'plate.twist.z',plateData.twist_z}; c=c+1;
            newData(c,:) = {'plate.twist.threeD',plateData.twist_3D}; c=c+1;
            newData(c,:) = {'plate.sisterCentreSpeed',plateData.sisCentreSpeed}; c=c+1;

            % get plate thickness measurements
            sisCentre_x = [];
            if length(dS.sisterList) == 1
              plateThickness = nan(1,nFrames);
            else
              for jSis = 1:length(dS.sisterList)
                sisCentre_x = [sisCentre_x nanmean([dS.sisterList(jSis).coords1(:,1) dS.sisterList(jSis).coords2(:,1)],2)];
              end
              plateThickness = nanstd(sisCentre_x,[],2);
            end
            % put data into string format
            newData(c,:) = {'plate.plateThickness',plateThickness}; c=c+1;

          
          % get intensities if required
          intsMean = nan(nFrames,2);
          intsMax = nan(nFrames,2);
          if isfield(dS,'spotInt')
              if size(sI.intensity,2)==1
                  for iFrame = 1:nFrames
                    % inner mean int
                    temp = cat(2,sI(iFrame).intensity) - bg(iFrame,1);
                    intsMean(iFrame,~isnan(spotIDs)) = temp(spotIDs(~isnan(spotIDs)),1);
                    % inner max int
                    temp = cat(2,sI(iFrame).intensity_max) - bg(iFrame,1);
                    intsMax(iFrame,~isnan(spotIDs)) = temp(spotIDs(~isnan(spotIDs)),1);
                  end
              else
                  for iFrame = 1:nFrames
                    % inner mean int
                    temp = cat(2,sI(iFrame).intensity) - bg(iFrame,1);
                    intsMean(iFrame,~isnan(spotIDs)) = temp(spotIDs(~isnan(spotIDs)),ch);
                    % inner max int
                    temp = cat(2,sI(iFrame).intensity_max) - bg(iFrame,1);
                    intsMax(iFrame,~isnan(spotIDs)) = temp(spotIDs(~isnan(spotIDs)),ch);
                  end
              end
          end %ints
          
          % put data into string format
          newData(c,:) = {'intensity.mean',intsMean}; c=c+1;
          newData(c,:) = {'intensity.max',intsMax}; c=c+1;
          newData(c,:) = {'intensity.bg',repmat(bg,nFrames,1)}; c=c+1;
            
            %% Directional information
            
            % get direction of movement
            direc = [];
            for iTrack = 1:2
              if ~isnan(trackIDs(iTrack)) && ~isempty(dS.trackList(trackIDs(iTrack)).direction)
                direc(:,iTrack) = dS.trackList(trackIDs(iTrack)).direction;
              else
                direc(:,iTrack) = nan(nFrames,1);
              end
            end
            direc(end+1,:) = NaN;
            
            % get potential switch events (i.e. individual timepoints between P and AP)
            switchBuffer = 4;
            switchEvent = zeros(size(direc));
            switchDirec = [diff(direc); NaN NaN];
            for iPoint = 1:2:size(switchEvent,1)-(switchBuffer+1)
                for jSis = 1:2
                    tempDirec = switchDirec(iPoint:iPoint+(switchBuffer-1),jSis);
                    idx = find(tempDirec(2:(switchBuffer-1))==0);
                    if max(abs(tempDirec))==1 && abs(sum(tempDirec))>1 && ~isempty(idx)
                        switchEvent(iPoint+idx(1):iPoint+idx(end),jSis) = 1;
                    end
                end
            end
            % calculate directional information
            direc_P = +(direc==1);               direc_P(direc_P==0) = NaN;
            direc_AP = +(direc==-1);             direc_AP(direc_AP==0) = NaN;
            direc_S = +(switchEvent==1);         direc_S(direc_S==0) = NaN;
            direc_N = +((direc+switchEvent)==0); direc_N(direc_N==0) = NaN;
            % put data into string format
            dirLabel = {'P','AP','S','N'};
            for iDir = 1:4
                eval(['newData(c,:) = {''direction.' dirLabel{iDir} ''',direc_' dirLabel{iDir} '}; c=c+1;']);
            end
        
            % compile new data with original
            allData = combineStrForms(allData,newData);
            
            % clear some data to ensure no overlap on next loop
            clear spotIDs newData direc
        
        end % sisters
    end % movies     
end % expts

%% Save results to structure

compiledInter = strForm2struct(allData);

%% Output any error information

if ~isempty(noSkip)
  fprintf('\nThe following cells were skipped by the user:\n');
  for iCell = 1:size(noSkip,1)
    fprintf('    Exp %i, Mov %i\n',noSkip(iCell,1),noSkip(iCell,2));
  end
end
if ~isempty(noFail)
  fprintf('\nThe following cells failed during spot detection:\n');
  for iCell = 1:size(noFail,1)
    fprintf('    Exp %i, Mov %i\n',noFail(iCell,1),noFail(iCell,2));
  end
end
if ~isempty(noFail)
  fprintf('\nThe following cells found no spots:\n');
  for iCell = 1:size(noFail,1)
    fprintf('    Exp %i, Mov %i\n',noFail(iCell,1),noFail(iCell,2));
  end
end
if ~isempty(noSis)
  fprintf('\nThe following cells contain no sisterList:\n');
  for iCell = 1:size(noSis,1)
    fprintf('    Exp %i, Mov %i\n',noSis(iCell,1),noSis(iCell,2));
  end
end
fprintf('\n');

end

function measurements = makeMeasurements(coords,plane)
    
    if nargin<2 || isempty(plane)
        plane = 0;
    end
        
    %% Inter-kinetochore: separation, twist, and velocities
    
    % coordinate-specific sister separation (using inner-kinetochore)
    sisSep_xyz = diff(reshape(coords,size(coords,1),3,2),[],3);
    measurements.sisSep_x = sisSep_xyz(:,1);
    measurements.sisSep_y = sisSep_xyz(:,2);
    measurements.sisSep_z = sisSep_xyz(:,3);

    % 3D sister separation
    sisSep_3D = sqrt(sum(sisSep_xyz.^2,2));
    measurements.sisSep_3D = sisSep_3D;

    % 2D sister separation
    sisSep_2D = sqrt(sum(sisSep_xyz(:,1:2).^2,2));
    measurements.sisSep_2D = sisSep_2D;

    if plane

        % coordinate-specific twist
        twist_y = sisSep_xyz(:,2)./sisSep_xyz(:,1);
        twist_y = atand(twist_y);
        twist_z = sisSep_xyz(:,3)./sisSep_xyz(:,1);
        twist_z = atand(twist_z);
        measurements.twist_y = twist_y;
        measurements.twist_z = twist_z;

        % 3D twist (dot product with the x-axis with length 1)
        xAxis = repmat([1 0 0],size(sisSep_xyz,1),1);
        twist_3D = dot(sisSep_xyz,xAxis,2);
        twist_3D = twist_3D./(sisSep_3D(:,1));
        twist_3D = acosd(twist_3D);
        twist_3D(twist_3D>90) = 180-twist_3D(twist_3D>90);
        measurements.twist_3D = twist_3D;

        % sister centre velocities in x-coordinate
        sisCentre_x = sum(coords(:,[1 4]),2)/2;
        sisCentreSpeed_x = [diff(sisCentre_x); NaN];
        measurements.sisCentreSpeed = sisCentreSpeed_x;

    end

end %makeMeasurements subfunction

