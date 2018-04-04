function compiledIntra = dublIntraMeasurements(movies,varargin)
% DUBLINTRAMEASUREMENTS Produces a structure of population-level
% intra-kinetochore measurements over multiple experiments.
%
%    DUBLINTRAMEASUREMENTS(MOVIES,...) Provides all coordinates, inter-
%    and intra-kinetochore measurements for all sisters within all cells
%    across all experiments. The resulting structure provides the tools for
%    deriving population-level analyses for an experiment. Options are
%    available.
%
%    Options, defaults in {}:-
%
%    centralise: {0} or 1. Whether or not to adjust the outer kinetochore
%       components' position so that average delta-x, delta-y and delta-z
%       measurements are zero.
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. The channels between
%       which to make intra-kinetochore measurements. The direction of
%       measurements will be defined by the channel orientation in the
%       neighbourSpots section of options.
%
%    intRefMarker: {'self'}, 'inner' or 'outer'. The marker around which
%       intensity measurements are measured. 'self' means intensity
%       measurements are measured about it's own spot centre, while 'inner'
%       and 'outer' mean intensity measurements are measured about the
%       inner or outer kinetochore marker, respectively.
%
%    paired: 0 or {1}. Whether or not to take paired measurements, or raw
%       spot-by-spot measurements.
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
opts.centralise = 1;
opts.channels = [1 2];
opts.intRefMarker = 'self';
opts.paired = 1;
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
if opts.paired
  for iExpt = 1:numExpts
    for iMov = 1:length(movies{iExpt})
      paired = isfield(movies{iExpt}{iMov}.dataStruct{opts.channels(1)},'sisterList');
      if paired; break; end
    end
    if paired; break; end
  end
end
if selType==3
    if paired
        subset = opts.spotSelection.selection;
        selType = 1;
    else
        subset = opts.spotSelection.rawSelection;
    end
end

if isempty(opts.prevMeas)
    
    % make new intra-measurements structure if no previous measurements
    % provided
    allData = dublMakeIntraStructure(paired);
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
    
    % find channel vector
    chanVect = movies{iExpt}{1}.options.neighbourSpots.channelOrientation;
    remChans = setdiff(chanVect,opts.channels);
    chanVect(chanVect==remChans) = [];
    
    for iMov = 1:length(theseMovies)
        
      % get the movie index
      movNum = theseMovies{iMov}.index;
      % check whether there is data in this movie
      if ~isfield(theseMovies{iMov},'dataStruct')
        noFail = [noFail; iExpt movNum];
        continue
      end
      
      % get dataStructs
      dSinner = theseMovies{iMov}.dataStruct{chanVect(1)};
      dSouter = theseMovies{iMov}.dataStruct{chanVect(2)};
      
      % check whether the movie failed
      if (isfield(dSinner,'failed') && dSinner.failed) || (isfield(dSouter,'failed') && dSouter.failed) ...
              || (~isfield(dSinner,'initCoord') && ~isfield(dSouter,'initCoord'))
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
      pixelSize = theseMovies{iMov}.metadata.pixelSize;
      % check if there is a plate fit
      plane = (isfield(dSinner,'planeFit') && ~isempty(dSinner.planeFit) && ~isempty(dSinner.planeFit.planeVectors));
      
      % get initCoord structures - if one doesn't exist trade for an empty
      % initCoord
      if ~isfield(dSinner,'initCoord')
        iCouter = dSouter.initCoord;
        iC = iCouter;
        iC.allCoord(:) = NaN; iC.allCoordPix(:) = NaN; iC.amp(:) = NaN; iC.bg(:) = NaN;
        iCinner = iC;
        clear iC
      elseif ~isfield(dSouter,'initCoord')
        iCinner = dSinner.initCoord;
        iC = iCinner;
        iC.allCoord(:) = NaN; iC.allCoordPix(:) = NaN; iC.amp(:) = NaN; iC.bg(:) = NaN;
        iCouter = iC;
        clear iC
      else
        iCinner = dSinner.initCoord;  
        iCouter = dSouter.initCoord;
      end
      
      % get xyz correction for centralisation
      if opts.centralise
        deltaxyz = [];
        for iFrame = 1:nFrames
          deltaxyz = [deltaxyz; iCouter(iFrame).allCoord(:,1:3)-iCinner(iFrame).allCoord(:,1:3)];
        end
        mCentMeds = nanmedian(deltaxyz);
      else
        mCentMeds = [0 0 0];
      end
      
      % get spotInts
      ints = (isfield(dSinner,'spotInt') && isfield(dSouter,'spotInt'));
      if ints
        sIinner = dSinner.spotInt;
        innerBg = dSinner.cellInt.back;
        sIouter = dSouter.spotInt;
        outerBg = dSouter.cellInt.back;
      else
        innerBg = nan(nFrames,1);
        outerBg = nan(nFrames,1);
      end
      
      if paired
        
        % check whether a sisterList is present, and if it contains any sisters
        if ~isfield(dSinner,'sisterList') || isempty(dSinner.sisterList(1).trackPairs) || ...
                ~isfield(dSouter,'sisterList') || isempty(dSouter.sisterList(1).trackPairs)
          noSis = [noSis; iExpt movNum];
          continue
        end
        
        % if no sisters given, go through all sisters in movie
        switch selType
            case 0 % no spot selection
                iSubset = 1:length(theseMovies{iMov}.dataStruct{chanVect(1)}.sisterList);
            case 1 % using spots/tracks
                iSubset = 1:length(theseMovies{iMov}.dataStruct{chanVect(1)}.sisterList);
                theseTracks = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
            case 2 % using sisters
                iSubset = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
        end

        % check that there are sisters
        if isempty(dSinner.sisterList(1).trackPairs)
            continue
        end
        
        for iSis = iSubset
            
            % construct sister pair label
            label = sprintf('%02d%02d%03d',iExpt,iMov,iSis);
            
            % start counter for storing data
            c=1;

            % get sisterLists
            sLinner = dSinner.sisterList(iSis);
            sLouter = dSouter.sisterList(iSis);
            
            % get trackID and spotIDs, make spotIDs nan if deselected
            trackIDs = dSinner.sisterList(1).trackPairs(iSis,1:2);
            spotIDs = nan(nFrames,2);
            for iTrack = 1:2
                if selType~=1 %none or sisters
                    spotIDs(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).featIndx;
                else %tracks
                    if ismember(trackIDs(iTrack),theseTracks)
                        spotIDs(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).featIndx;
                    else
                        trackIDs(iTrack) = NaN;
                        switch iTrack
                            case 1
                                sLinner.coords1(:,1:3) = NaN;
                                sLouter.coords1(:,1:3) = NaN;
                            case 2
                                sLinner.coords2(:,1:3) = NaN;
                                sLouter.coords2(:,1:3) = NaN;
                        end
                    end
                end 
            end
            % if both spots skipped
            if all(isnan(trackIDs))
                continue;
            end

            %% Kinetochore positions
            
            % get microscope coordinates, and plane coordinates if present
            mCoordsInner = nan(nFrames,6);
            mCoordsOuter = nan(nFrames,6);
            pCoordsInner = nan(nFrames,6);
            pCoordsOuter = nan(nFrames,6);
            for iFrame = 1:nFrames
                for iTrk = 1:2
                    % track one stored in (:,1:3), two in (:,4:6)
                    rng = (3*(iTrk-1)+1):3*iTrk;
                    if ~isnan(spotIDs(iFrame,iTrk))
                        mCoordsInner(iFrame,rng) = iCinner(iFrame).allCoord(spotIDs(iFrame,iTrk),1:3);
                        mCoordsOuter(iFrame,rng) = iCouter(iFrame).allCoord(spotIDs(iFrame,iTrk),1:3) - mCentMeds;
                        % check whether or not this movie has a planeFit
                        if plane
                            % rotate coordinates into plane
                            pF = dSinner.planeFit;
                            if ~isempty(pF(iFrame).planeVectors)
                                coordSystem = pF(iFrame).planeVectors;
                                pCoordsInner(iFrame,:) = mCoordsInner(iFrame,:) - repmat(pF.planeOrigin,1,2);
                                pCoordsOuter(iFrame,:) = mCoordsOuter(iFrame,:) - repmat(pF.planeOrigin,1,2);
                                
                                pCoordsInner(iFrame,1:3) = (coordSystem\(pCoordsInner(iFrame,1:3)'))';
                                pCoordsInner(iFrame,4:6) = (coordSystem\(pCoordsInner(iFrame,4:6)'))';
                                pCoordsOuter(iFrame,1:3) = (coordSystem\(pCoordsOuter(iFrame,1:3)'))';
                                pCoordsOuter(iFrame,4:6) = (coordSystem\(pCoordsOuter(iFrame,4:6)'))';
                            end
                        end
                    end
                end
            end
            
            % put data into string format
            newData(c,:) = {'label',label}; c=c+1;
            
            % get microscope coordinates of each spot
            mCoords_x = [mCoordsInner(:,[1 4]) mCoordsOuter(:,[1 4])];
            mCoords_y = [mCoordsInner(:,[2 5]) mCoordsOuter(:,[2 5])];
            mCoords_z = [mCoordsInner(:,[3 6]) mCoordsOuter(:,[3 6])];
            % put data into string format
            newData(c,:) = {'microscope.coords.x',mCoords_x(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'microscope.coords.y',mCoords_y(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'microscope.coords.z',mCoords_z(:,[1 3 2 4])}; c=c+1;
            
            % get plate coordinates of each spot
            pCoords_x = [pCoordsInner(:,[1 4]) pCoordsOuter(:,[1 4])];
            pCoords_y = [pCoordsInner(:,[2 5]) pCoordsOuter(:,[2 5])];
            pCoords_z = [pCoordsInner(:,[3 6]) pCoordsOuter(:,[3 6])];
            % put data into string format
            newData(c,:) = {'plate.coords.x',pCoords_x(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'plate.coords.y',pCoords_y(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'plate.coords.z',pCoords_z(:,[1 3 2 4])}; c=c+1;
            
            %% Inter- and intra-kinetochore measurements
            
            micrData = pairedMeasurements(mCoordsInner,mCoordsOuter,0);
            % put data into string format
            newData(c,:) = {'microscope.sisSep.x',micrData.sisSep_x}; c=c+1;
            newData(c,:) = {'microscope.sisSep.y',micrData.sisSep_y}; c=c+1;
            newData(c,:) = {'microscope.sisSep.z',micrData.sisSep_z}; c=c+1;
            newData(c,:) = {'microscope.sisSep.twoD',micrData.sisSep_2D}; c=c+1;
            newData(c,:) = {'microscope.sisSep.threeD',micrData.sisSep_3D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.x.all',micrData.delta_x}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.y.all',micrData.delta_y}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.z.all',micrData.delta_z}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.oneD',micrData.delta_1D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.twoD.all',micrData.delta_2D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.threeD.all',micrData.delta_3D}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.y.all',micrData.swivel_y}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.z.all',micrData.swivel_z}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.threeD.all',micrData.swivel_3D}; c=c+1;
            if isfield(micrData,'swivel_kMT')
              newData(c,:) = {'microscope.raw.swivel.kMT',micrData.swivel_kMT}; c=c+1;
            else
              newData(c,:) = {'microscope.raw.swivel.kMT',[]}; c=c+1;
            end
            
            plateData = pairedMeasurements(pCoordsInner,pCoordsOuter,1);
            % put data into string format
            newData(c,:) = {'plate.sisSep.x',plateData.sisSep_x}; c=c+1;
            newData(c,:) = {'plate.sisSep.y',plateData.sisSep_y}; c=c+1;
            newData(c,:) = {'plate.sisSep.z',plateData.sisSep_z}; c=c+1;
            newData(c,:) = {'plate.sisSep.twoD',plateData.sisSep_2D}; c=c+1;
            newData(c,:) = {'plate.sisSep.threeD',plateData.sisSep_3D}; c=c+1;
            newData(c,:) = {'plate.raw.delta.x.all',plateData.delta_x}; c=c+1;
            newData(c,:) = {'plate.raw.delta.y.all',plateData.delta_y}; c=c+1;
            newData(c,:) = {'plate.raw.delta.z.all',plateData.delta_z}; c=c+1;
            newData(c,:) = {'plate.raw.delta.oneD',plateData.delta_1D}; c=c+1;
            newData(c,:) = {'plate.raw.delta.twoD.all',plateData.delta_2D}; c=c+1;
            newData(c,:) = {'plate.raw.delta.threeD.all',plateData.delta_3D}; c=c+1;
            newData(c,:) = {'plate.raw.swivel.y.all',plateData.swivel_y}; c=c+1;
            newData(c,:) = {'plate.raw.swivel.z.all',plateData.swivel_z}; c=c+1;
            newData(c,:) = {'plate.raw.swivel.threeD.all',plateData.swivel_3D}; c=c+1;
            if isfield(plateData,'swivel_kMT')
              newData(c,:) = {'plate.raw.swivel.kMT',plateData.swivel_kMT}; c=c+1;
            else
              newData(c,:) = {'plate.raw.swivel.kMT',[]}; c=c+1;
            end
            newData(c,:) = {'plate.twist.y',plateData.twist_y}; c=c+1;
            newData(c,:) = {'plate.twist.z',plateData.twist_z}; c=c+1;
            newData(c,:) = {'plate.twist.threeD',plateData.twist_3D}; c=c+1;
            newData(c,:) = {'plate.sisterCentreSpeed',plateData.sisCentreSpeed}; c=c+1;

            % get plate thickness measurements
            sisCentre_x = [];
            if length(dSinner.sisterList) == 1
              plateThickness = nan(1,nFrames);
            else
              for jSis = 1:length(dSinner.sisterList)
                sisCentre_x = [sisCentre_x nanmean([dSinner.sisterList(jSis).coords1(:,1) dSinner.sisterList(jSis).coords2(:,1)],2)];
              end
              plateThickness = nanstd(sisCentre_x,[],2);
            end
            % put data into string format
            newData(c,:) = {'plate.plateThickness',plateThickness}; c=c+1;

            %% Quality control
            
            % find which data satisfies the z-depth requirement
            satisfies = +(abs(micrData.delta_z)<0.5*pixelSize(3));
            satisfies(satisfies==0) = NaN;
            % put all data into string format
            newData(c,:) = {'microscope.depthFilter.delta.x.all',micrData.delta_x.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.y.all',micrData.delta_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.z.all',micrData.delta_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.oneD',micrData.delta_1D.*prod(satisfies,2)}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.twoD.all',micrData.delta_2D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.threeD.all',micrData.delta_3D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.y.all',micrData.swivel_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.z.all',micrData.swivel_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.threeD.all',micrData.swivel_3D.*satisfies}; c=c+1;
            if isfield(micrData,'swivel_kMT')
              newData(c,:) = {'microscope.depthFilter.swivel.kMT',micrData.swivel_kMT.*satisfies}; c=c+1;
            else
              newData(c,:) = {'microscope.depthFilter.swivel.kMT',[]}; c=c+1;
            end
            newData(c,:) = {'plate.depthFilter.delta.x.all',plateData.delta_x.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.delta.y.all',plateData.delta_y.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.delta.z.all',plateData.delta_z.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.delta.oneD',plateData.delta_1D.*prod(satisfies,2)}; c=c+1;
            newData(c,:) = {'plate.depthFilter.delta.twoD.all',plateData.delta_2D.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.delta.threeD.all',plateData.delta_3D.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.swivel.y.all',plateData.swivel_y.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.swivel.z.all',plateData.swivel_z.*satisfies}; c=c+1;
            newData(c,:) = {'plate.depthFilter.swivel.threeD.all',plateData.swivel_3D.*satisfies}; c=c+1;
            if isfield(plateData,'swivel_kMT')
              newData(c,:) = {'plate.depthFilter.swivel.kMT',plateData.swivel_kMT.*satisfies}; c=c+1;
            else
              newData(c,:) = {'plate.depthFilter.swivel.kMT',[]}; c=c+1;
            end
            
          % get intensities if required
          intsInnerMean = nan(nFrames,2);
          intsOuterMean = nan(nFrames,2);
          intsInnerMax = nan(nFrames,2);
          intsOuterMax = nan(nFrames,2);
          if ints
            if size(sIinner.intensity,2)==1
              for iFrame = 1:nFrames
                % inner mean int
                temp = cat(2,sIinner(iFrame).intensity) - innerBg(iFrame,1);
                intsInnerMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),1);
                % inner max int
                temp = cat(2,sIinner(iFrame).intensity_max) - innerBg(iFrame,1);
                intsInnerMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),1);
                % outer mean int
                temp = cat(2,sIouter(iFrame).intensity) - outerBg(iFrame,1);
                intsOuterMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),1);
                % outer max int
                temp = cat(2,sIouter(iFrame).intensity_max) - outerBg(iFrame,1);
                intsOuterMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),1);
              end
            else
              switch opts.intRefMarker
                case 'self'
                  for iFrame = 1:nFrames
                    % inner mean int
                    temp = cat(2,sIinner(iFrame).intensity) - innerBg(iFrame,1);
                    intsInnerMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(1));
                    % inner max int
                    temp = cat(2,sIinner(iFrame).intensity_max) - innerBg(iFrame,1);
                    intsInnerMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(1));
                    % outer mean int
                    temp = cat(2,sIouter(iFrame).intensity) - outerBg(iFrame,1);
                    intsOuterMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(2));
                    % outer max int
                    temp = cat(2,sIouter(iFrame).intensity_max) - outerBg(iFrame,1);
                    intsOuterMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(2));
                  end
                case 'inner'
                  for iFrame = 1:nFrames
                    % inner mean int
                    temp = cat(2,sIinner(iFrame).intensity) - innerBg(iFrame,1);
                    intsInnerMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(1));
                    % inner max int
                    temp = cat(2,sIinner(iFrame).intensity_max) - innerBg(iFrame,1);
                    intsInnerMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(1));
                    % outer mean int, relative to inner
                    temp = cat(2,sIinner(iFrame).intensity) - outerBg(iFrame,1);
                    intsOuterMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(2));
                    % outer max int, relative to inner
                    temp = cat(2,sIinner(iFrame).intensity_max) - outerBg(iFrame,1);
                    intsOuterMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(2));
                  end
                case 'outer'
                  for iFrame = 1:nFrames
                    % inner mean int, relative to outer
                    temp = cat(2,sIouter(iFrame).intensity) - innerBg(iFrame,1);
                    intsInnerMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(1));
                    % inner max int, relative to outer
                    temp = cat(2,sIouter(iFrame).intensity_max) - innerBg(iFrame,1);
                    intsInnerMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(1));
                    % outer mean int
                    temp = cat(2,sIouter(iFrame).intensity) - outerBg(iFrame,1);
                    intsOuterMean(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(2));
                    % outer max int
                    temp = cat(2,sIouter(iFrame).intensity_max) - outerBg(iFrame,1);
                    intsOuterMax(iFrame,:) = temp(spotIDs(~isnan(spotIDs)),chanVect(2));
                  end
              end
            end
          end %ints
          
          % put data into string format
          newData(c,:) = {'intensity.mean.inner',intsInnerMean}; c=c+1;
          newData(c,:) = {'intensity.mean.outer',intsOuterMean}; c=c+1;
          newData(c,:) = {'intensity.max.inner',intsInnerMax}; c=c+1;
          newData(c,:) = {'intensity.max.outer',intsOuterMax}; c=c+1;
          newData(c,:) = {'intensity.bg.inner',repmat(innerBg,nFrames,1)}; c=c+1;
          newData(c,:) = {'intensity.bg.outer',repmat(outerBg,nFrames,1)}; c=c+1;
            
            %% Directional information
            
            % get direction of movement
            direc = [];
            for iTrack = 1:2
              if ~isnan(trackIDs(iTrack)) && ~isempty(dSinner.trackList(trackIDs(iTrack)).direction)
                direc(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).direction;
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
        
      else
          
          % start counter for storing data
          c=1;
          
          % if no sisters given, go through all sisters in movie
          switch selType
            case 0 % no spot selection
                spotIDs = 1:size(iCinner(1).allCoord,1);
            case 1 % using spots/tracks
                trackIDs = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
                spotIDs = cat(2,dSinner.trackList(trackIDs).featIndx);
            case 3 % using initCoord
                spotIDs = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
          end
          nSpots = length(spotIDs);
          if nSpots == 0
              continue
          end
          
          % construct spot label
          for iSpot = 1:nSpots
            labels(iSpot,:) = sprintf('%02d%02d%03d',iExpt,iMov,iSpot);
          end
          
          % put data into string format
          newData(c,:) = {'label',labels}; c=c+1;
          
          %% Kinetochore positions
          
          % get microscope coordinates, and plane coordinates if present
          mCoordsInner = iCinner(1).allCoord(spotIDs,1:3);
          mCoordsOuter = iCouter(1).allCoord(spotIDs,1:3) - repmat(mCentMeds,nSpots,1);
          % check whether or not this movie has a planeFit
          if plane
              % rotate coordinates into plane
              pF = dSinner.planeFit;
              if ~isempty(pF(1).planeVectors)
                  coordSystem = pF(1).planeVectors;
                  pCoordsInner = (coordSystem\(mCoordsInner'))';
                  pCoordsOuter = (coordSystem\(mCoordsOuter'))';
              end
          else
              % give empty datasets the size of microscopy coordinates
              pCoordsInner = nan(size(mCoordsInner));
              pCoordsOuter = nan(size(mCoordsOuter));
          end
            
          % get microscope coordinates of each spot
          mCoords_x = [mCoordsInner(:,1) mCoordsOuter(:,1)];
          mCoords_y = [mCoordsInner(:,2) mCoordsOuter(:,2)];
          mCoords_z = [mCoordsInner(:,3) mCoordsOuter(:,3)];
          % put data into string format
          newData(c,:) = {'microscope.coords.x',mCoords_x}; c=c+1;
          newData(c,:) = {'microscope.coords.y',mCoords_y}; c=c+1;
          newData(c,:) = {'microscope.coords.z',mCoords_z}; c=c+1;
          
          % get plate coordinates of each spot
          pCoords_x = [pCoordsInner(:,1) pCoordsOuter(:,1)];
          pCoords_y = [pCoordsInner(:,2) pCoordsOuter(:,2)];
          pCoords_z = [pCoordsInner(:,3) pCoordsOuter(:,3)];
          % put data into string format
          newData(c,:) = {'plate.coords.x',pCoords_x}; c=c+1;
          newData(c,:) = {'plate.coords.y',pCoords_y}; c=c+1;
          newData(c,:) = {'plate.coords.z',pCoords_z}; c=c+1;
          
          % get intensities if required
          intsInnerMean = nan(nSpots,1);
          intsOuterMean = nan(nSpots,1);
          intsInnerMax = nan(nSpots,1);
          intsOuterMax = nan(nSpots,1);
          if ints
            if size(sIinner.intensity,2)==1
              % inner mean int
              temp = cat(2,sIinner.intensity) - innerBg;
              intsInnerMean = temp(spotIDs,1);
              % inner max int
              temp = cat(2,sIinner.intensity_max) - innerBg;
              intsInnerMax = temp(spotIDs,1);
              % outer mean int
              temp = cat(2,sIouter.intensity) - outerBg;
              intsOuterMean = temp(spotIDs,1);
              % outer max int
              temp = cat(2,sIouter.intensity_max) - outerBg;
              intsOuterMax = temp(spotIDs,1);
            else
              switch opts.intRefMarker
                case 'self'
                  % inner mean int
                  temp = cat(2,sIinner.intensity) - innerBg;
                  intsInnerMean = temp(spotIDs,chanVect(1));
                  % inner max int
                  temp = cat(2,sIinner.intensity_max) - innerBg;
                  intsInnerMax = temp(spotIDs,chanVect(1));
                  % outer mean int
                  temp = cat(2,sIouter.intensity) - outerBg;
                  intsOuterMean = temp(spotIDs,chanVect(2));
                  % outer max int
                  temp = cat(2,sIouter.intensity_max) - outerBg;
                  intsOuterMax = temp(spotIDs,chanVect(2));
                case 'inner'
                  % inner mean int
                  temp = cat(2,sIinner.intensity) - innerBg;
                  intsInnerMean = temp(spotIDs,chanVect(1));
                  % inner max int
                  temp = cat(2,sIinner.intensity_max) - innerBg;
                  intsInnerMax = temp(spotIDs,chanVect(1));
                  % outer mean int, relative to inner
                  temp = cat(2,sIinner.intensity) - outerBg;
                  intsOuterMean = temp(spotIDs,chanVect(2));
                  % outer max int, relative to inner
                  temp = cat(2,sIinner.intensity_max) - outerBg;
                  intsOuterMax = temp(spotIDs,chanVect(2));
                case 'outer'
                  % inner mean int, relative to outer
                  temp = cat(2,sIouter.intensity) - innerBg;
                  intsInnerMean = temp(spotIDs,chanVect(1));
                  % inner max int, relative to outer
                  temp = cat(2,sIouter.intensity_max) - innerBg;
                  intsInnerMax = temp(spotIDs,chanVect(1));
                  % outer mean int
                  temp = cat(2,sIouter.intensity) - outerBg;
                  intsOuterMean = temp(spotIDs,chanVect(2));
                  % outer max int
                  temp = cat(2,sIouter.intensity_max) - outerBg;
                  intsOuterMax = temp(spotIDs,chanVect(2));
              end
            end
          end %ints
          
          % put data into string format
          newData(c,:) = {'intensity.mean.inner',intsInnerMean}; c=c+1;
          newData(c,:) = {'intensity.mean.outer',intsOuterMean}; c=c+1;
          newData(c,:) = {'intensity.max.inner',intsInnerMax}; c=c+1;
          newData(c,:) = {'intensity.max.outer',intsOuterMax}; c=c+1;
          newData(c,:) = {'intensity.bg.inner',repmat(innerBg,nSpots,1)}; c=c+1;
          newData(c,:) = {'intensity.bg.outer',repmat(outerBg,nSpots,1)}; c=c+1;
          
          %% Inter- and intra-kinetochore measurements
            
          micrData = soloMeasurements(mCoordsInner,mCoordsOuter);
          % put data into string format
          newData(c,:) = {'microscope.raw.delta.x.all',micrData.delta_x}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.y.all',micrData.delta_y}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.z.all',micrData.delta_z}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.twoD.all',micrData.delta_2D}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.threeD.all',micrData.delta_3D}; c=c+1;
          
          plateData = soloMeasurements(pCoordsInner,pCoordsOuter);
          % put data into string format
          newData(c,:) = {'plate.raw.delta.x.all',plateData.delta_x}; c=c+1;
          newData(c,:) = {'plate.raw.delta.y.all',plateData.delta_y}; c=c+1;
          newData(c,:) = {'plate.raw.delta.z.all',plateData.delta_z}; c=c+1;
          newData(c,:) = {'plate.raw.delta.twoD.all',plateData.delta_2D}; c=c+1;
          newData(c,:) = {'plate.raw.delta.threeD.all',plateData.delta_3D}; c=c+1;
          
          %% Quality control
            
          % find which data satisfies the z-depth requirement
          satisfies = +(abs(micrData.delta_z)<0.5*pixelSize(3));
          satisfies(satisfies==0) = NaN;
          
          % put all data into string format
          newData(c,:) = {'microscope.depthFilter.delta.x.all',micrData.delta_x.*satisfies}; c=c+1;
          newData(c,:) = {'microscope.depthFilter.delta.y.all',micrData.delta_y.*satisfies}; c=c+1;
          newData(c,:) = {'microscope.depthFilter.delta.z.all',micrData.delta_z.*satisfies}; c=c+1;
          newData(c,:) = {'microscope.depthFilter.delta.twoD.all',micrData.delta_2D.*satisfies}; c=c+1;
          newData(c,:) = {'microscope.depthFilter.delta.threeD.all',micrData.delta_3D.*satisfies}; c=c+1;
          
          newData(c,:) = {'plate.depthFilter.delta.x.all',plateData.delta_x.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.y.all',plateData.delta_y.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.z.all',plateData.delta_z.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.twoD.all',plateData.delta_2D.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.threeD.all',plateData.delta_3D.*satisfies}; c=c+1;
          
          % compile new data with original
          allData = combineStrForms(allData,newData);  
          clear newData spotIDs
            
      end % paired
    end % movies     
end % expts

%% Save results to structure

compiledIntra = strForm2struct(allData);

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

function measurements = pairedMeasurements(coordsInner,coordsOuter,plane)
    
    if nargin<3 || isempty(plane)
        plane = 0;
    end
        
    %% Inter-kinetochore: separation, twist, and velocities
    
    % coordinate-specific sister separation (using inner-kinetochore)
    sisSep_xyz = diff(reshape(coordsInner,size(coordsInner,1),3,2),[],3);
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
        sisCentre_x = sum(coordsInner(:,[1 4]),2)/2;
        sisCentreSpeed_x = [diff(sisCentre_x); NaN];
        measurements.sisCentreSpeed = sisCentreSpeed_x;

    end

    %% Intra-kinetochore delta vector
    
    % coordinate-specific delta
    delta1_xyz = coordsOuter(:,1:3) - coordsInner(:,1:3);
    delta2_xyz = coordsOuter(:,4:6) - coordsInner(:,4:6);
    
    measurements.delta_x = [delta1_xyz(:,1) delta2_xyz(:,1)];
    measurements.delta_y = [delta1_xyz(:,2) delta2_xyz(:,2)];
    measurements.delta_z = [delta1_xyz(:,3) delta2_xyz(:,3)];
    
    % 3D delta
    delta1_3D = sqrt(sum(delta1_xyz.^2,2));
    delta2_3D = sqrt(sum(delta2_xyz.^2,2));
    measurements.delta_3D = [delta1_3D delta2_3D];
    
    % 2D delta
    delta1_2D = sqrt(sum(delta1_xyz(:,1:2).^2,2));
    delta2_2D = sqrt(sum(delta2_xyz(:,1:2).^2,2));
    measurements.delta_2D = [delta1_2D delta2_2D];
    
    % 1D delta
    outerSisSep_xyz = diff(reshape(coordsOuter,size(coordsOuter,1),3,2),[],3);
    outerSisSep_2D = sqrt(sum(outerSisSep_xyz(:,1:2).^2,2));
    delta_1D = (outerSisSep_2D-sisSep_2D)/2;
    measurements.delta_1D = delta_1D;
    
    %% Intra-kinetochore swivel
    
    % y swivel
    dy = sisSep_2D;
    del = delta1_2D;
    eps = sqrt(sum((coordsInner(:,4:5)-coordsOuter(:,1:2)).^2,2));
    swivel1_y = (dy.^2+del.^2-eps.^2)./(2*dy.*del);
    swivel1_y = 180 - acosd(swivel1_y);
    swivel1_y = swivel1_y.*sign(delta1_xyz(:,2));

    del = delta2_2D;
    eps = sqrt(sum((coordsOuter(:,4:5)-coordsInner(:,1:2)).^2,2));
    swivel2_y = (dy.^2+del.^2-eps.^2)./(2*dy.*del);
    swivel2_y = 180 - acosd(swivel2_y);
    swivel2_y = swivel2_y.*sign(delta2_xyz(:,2));

    measurements.swivel_y = [swivel1_y swivel2_y];

    % z swivel
    dz = sqrt(sum(sisSep_xyz(:,[1 3]).^2,2));
    del = sqrt(sum((coordsInner(:,[1 3])-coordsOuter(:,[1 3])).^2,2));
    eps = sqrt(sum((coordsOuter(:,[1 3])-coordsInner(:,[4 6])).^2,2));
    swivel1_z = (dz.^2+del.^2-eps.^2)./(2*dz.*del);
    swivel1_z = 180 - acosd(swivel1_z);
    swivel1_z = swivel1_z.*sign(delta1_xyz(:,3));

    del = sqrt(sum((coordsInner(:,[4 6])-coordsOuter(:,[4 6])).^2,2));
    eps = sqrt(sum((coordsOuter(:,[4 6])-coordsInner(:,[1 3])).^2,2));
    swivel2_z = (dz.^2+del.^2-eps.^2)./(2*dz.*del);
    swivel2_z = 180 - acosd(swivel2_z);
    swivel2_z = swivel2_z.*sign(delta2_xyz(:,3));

    measurements.swivel_z = [swivel1_z swivel2_z];

    % 3D swivel
    swivel1_3D = dot(-sisSep_xyz,delta1_xyz,2);
    swivel1_3D = swivel1_3D./(sisSep_3D(:,1).*delta1_3D(:,1));
    swivel1_3D = acosd(swivel1_3D); 
    swivel2_3D = dot(sisSep_xyz,delta2_xyz,2);
    swivel2_3D = swivel2_3D./(sisSep_3D(:,1).*delta2_3D(:,1));
    swivel2_3D = acosd(swivel2_3D);
    measurements.swivel_3D = [swivel1_3D swivel2_3D];    

end %pairedMeasurements subfunction

function measurements = soloMeasurements(coordsInner,coordsOuter)
    
    % coordinate-specific delta
    delta_xyz = coordsOuter(:,1:3) - coordsInner(:,1:3);
    
    measurements.delta_x = delta_xyz(:,1);
    measurements.delta_y = delta_xyz(:,2);
    measurements.delta_z = delta_xyz(:,3);
    
    % 3D delta
    delta_3D = sqrt(sum(delta_xyz.^2,2));
    measurements.delta_3D = delta_3D;
    
    % 2D delta
    delta_2D = sqrt(sum(delta_xyz(:,1:2).^2,2));
    measurements.delta_2D = delta_2D;

end %soloMeasurements subfunction

end


