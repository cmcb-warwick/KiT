function compiledInts = dublIntensityMeasurements(movies,varargin)
% DUBLINTENSITYMEASUREMENTS Produces a structure of population-level
% kinetochore intensity measurements over multiple experiments.
%
%    DUBLINTENSITYMEASUREMENTS(MOVIES,...) Provides all kinetochore
%    intensity measurements for all sisters within all cells across all
%    experiments. The resulting structure provides the tools for deriving
%    population-level intensity analyses for an experiment. Options are 
%    available.
%
%    Options, defaults in {}:-
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. The channels between
%       which to make intensity measurements.
%
%    refMarker: {'self'}, 'inner' or 'outer'. The marker around which
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
opts.channels = [1 2];
opts.refMarker = 'self';
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

% find whether any movies have paired spots
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
    allData = dublMakeIntensityStructure();
    allData = struct2strForm(allData);
    
else
    % get all old data
    allData = struct2strForm(opts.prevMeas);
    
end

% predesignate error arrays
noDS = [];
noSpot = [];
noSis = [];

%% Compiling measurements

for iExpt = 1:numExpts
    
    % get movie and sister list
    theseMovies = movies{iExpt};
    iSubset = subset{iExpt};
    
    % find channel vector
    chanVect = movies{iExpt}{1}.options.neighbourSpots.channelOrientation;
    chanVect = intersect(chanVect,opts.channels,'stable');
    
    % find spot reference channel
    refChan = movies{iExpt}{1}.options.coordSystemChannel;
    
    for iMov = 1:length(theseMovies)
      
      % get the movie index
      movNum = theseMovies{iMov}.index;
      % check whether there is data in this movie
      if ~isfield(theseMovies{iMov},'dataStruct')
        noDS = [noDS; iExpt movNum];
        continue
      end
      
      % get dataStructs
      dSinner = theseMovies{iMov}.dataStruct{chanVect(1)};
      dSouter = theseMovies{iMov}.dataStruct{chanVect(2)};
      refdS   = theseMovies{iMov}.dataStruct{refChan};
      
      % check whether the movie failed
      if ~isfield(refdS,'failed') || refdS.failed || ~isfield(refdS,'failed') || refdS.failed
        noSpot = [noSpot; iExpt movNum];
        continue
      end
      
      % get basic metadata
      nFrames = theseMovies{iMov}.metadata.nFrames;
      
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
        if ~isfield(refdS,'sisterList') || isempty(refdS.sisterList(1).trackPairs) || ...
                ~isfield(refdS,'sisterList') || isempty(refdS.sisterList(1).trackPairs)
          noSis = [noSis; iExpt movNum];
          continue
        end
        
        % if no sisters given, go through all sisters in movie
        switch selType
            case 0 % no spot selection
                iSubset = 1:length(refdS.sisterList);
            case 1 % using spots/tracks
                iSubset = 1:length(refdS.sisterList);
                theseTracks = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
            case 2 % using sisters
                iSubset = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
        end
        
        for iSis = iSubset
            
          % construct sister pair label
          label = sprintf('%02d%02d%03d',iExpt,iMov,iSis);
            
          % start counter for storing data
          c=1;

          % get trackID and spotIDs, make spotIDs nan if deselected
          trackIDs = refdS.sisterList(1).trackPairs(iSis,1:2);
          spotIDs = nan(nFrames,2);
          for iTrack = 1:2
            if selType~=1 %none or sisters
              spotIDs(:,iTrack) = refdS.trackList(trackIDs(iTrack)).featIndx;
            else %tracks
              if ismember(trackIDs(iTrack),theseTracks)
                spotIDs(:,iTrack) = refdS.trackList(trackIDs(iTrack)).featIndx;
              else
                trackIDs(iTrack) = NaN;
              end
            end 
          end
          % if both spots skipped
          if all(isnan(trackIDs))
            continue
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
              switch opts.refMarker
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
          newData(c,:) = {'intensity.bg.inner',repmat(innerBg,size(intsInnerMean,1),1)}; c=c+1;
          newData(c,:) = {'intensity.bg.outer',repmat(outerBg,size(intsOuterMean,1),1)}; c=c+1;
        
          % compile new data with original
          allData = combineStrForms(allData,newData);
          
          % clear some data to ensure no overlap on next loop
          clear spotIDs newData
        
        end % sisters
        
      else
          
          % start counter for storing data
          c=1;
          
          % if no sisters given, go through all sisters in movie
          switch selType
            case 0 % no spot selection
                spotIDs = 1:size(refdS.initCoord(1).allCoord,1);
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
              switch opts.refMarker
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
          newData(c,:) = {'intensity.bg.inner',repmat(innerBg,size(intsInnerMean,1),1)}; c=c+1;
          newData(c,:) = {'intensity.bg.outer',repmat(outerBg,size(intsOuterMean,1),1)}; c=c+1;
          
          % compile new data with original
          allData = combineStrForms(allData,newData);  
          clear newData spotIDs
            
      end % paired
    end % movies     
end % expts

%% Save results to structure

compiledInts = strForm2struct(allData);

%% Output any error information

if ~isempty(noDS)
  fprintf('\nThe following cells failed during spot detection:\n');
  for iCell = 1:size(noDS,1)
    fprintf('    Exp %i, Mov %i\n',noDS(iCell,1),noDS(iCell,2));
  end
end
if ~isempty(noSpot)
  fprintf('\nThe following cells found no spots:\n');
  for iCell = 1:size(noSpot,1)
    fprintf('    Exp %i, Mov %i\n',noSpot(iCell,1),noSpot(iCell,2));
  end
end
if ~isempty(noSis)
  fprintf('\nThe following cells contain no sisterList:\n');
  for iCell = 1:size(noSis,1)
    fprintf('    Exp %i, Mov %i\n',noSis(iCell,1),noSis(iCell,2));
  end
end
fprintf('\n');

