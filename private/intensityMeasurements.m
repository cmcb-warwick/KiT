function varargout = intensityMeasurements(movieStructs,varargin)
% Produce an analysis structure pertaining to intensity measurements.
%
%

% default vs user-defined options
opts.analysisChans = [1 2];
opts.controlChan = 2;
opts.measure = 'mean'; % mean, max, median, min
opts = processOptions(opts,varargin{:});

% process input
if nargin<1 || isempty(movieStructs)
    error('Please provide a movie structure.')
elseif ~iscell(movieStructs{1})
    movieStructs = {movieStructs};
end
if length(opts.controlChan)>1
    opts.controlChan = opts.controlChan(1);
    kitLog('Cannot have multiple control channels. Using channel %i as control.',opts.controlChan);
end
if ismember(opts.controlChan,opts.analysisChans)
    % remove control chan from list of analysis channels
    opts.channels = opts.analysisChans;
    opts.analysisChans = setdiff(opts.analysisChans,opts.controlChan);
else
    opts.channels = sort([opts.controlChan opts.analysisChans]);
end
opts.analysisChans = opts.analysisChans(:)';

% get some basic information
nExpts = length(movieStructs);
norms = struct('cellwise',[],'spotwise',[]);
analysis = struct('raw',[],...
                  'norm',norms,...
                  'bg',[],...
                  'stats',[]);

%% Collect data from movies

% make collection structures
allSpotInts = [];
allBg = [];
allNorms = [];
for iExpt = 1:nExpts
    
    % get this movieStruct
    thisMS = movieStructs{iExpt};
    nMovs = length(thisMS);
    
    for iMov = 1:nMovs
        
        % get this movie
        thisMov = thisMS{iMov};
        
        % check whether there is any data in the given channels
        if ~isfield(thisMov,'dataStruct')
            kitLog('No dataStruct detected: expt %i, mov %i.', iExpt, iMov);
            continue
        elseif length(thisMov.dataStruct) < max(opts.channels)
            kitLog('Data structure not present for all analysis channels: expt %i, mov %i.', iExpt, iMov);
            continue
        elseif any(~ismember(opts.channels,find(thisMov.options.intensity.execute == 1)))
            kitLog('Not all analysis channels were measured for intensities: expt %i, mov %i.', iExpt, iMov);
            continue
        elseif ~isfield(thisMov.dataStruct{opts.controlChan},'spotInt')
            kitLog('No spots found from which to measure intensities: expt %i, mov %i.', iExpt, iMov);
            continue
        end
        
        % set up collection arrays
        nSpots = length(thisMov.dataStruct{opts.controlChan}.spotInt);
        nTPs = length(thisMov.dataStruct{opts.controlChan}.spotInt(1).intensity);
        spotInts = nan(nSpots*nTPs,3);
        bg = nan(nTPs,3);
        
        % loop over channels
        for iChan = opts.channels
            
            % get this channel's spot intensities
            thisSpotInt = thisMov.dataStruct{iChan}.spotInt;
            thisCellInt = thisMov.dataStruct{iChan}.cellInt;
            % get background intensities
            bg(:,iChan) = thisCellInt.back(:);
            
            % THIS SHOULD BE REMOVED LATER - if number of spots between
            % channels differ, skip
            if length(thisSpotInt) ~= nSpots
              continue
            end
            
            % loop over spots to get their intensities
            for iSpot = 1:nSpots
                
                % check whether or not an intensity has been measured
                if isempty(thisSpotInt) || any(thisSpotInt(iSpot).intensity==[0,NaN])
                    continue
                end
                
                % find print range
                rnge = nTPs*(iSpot-1)+1:nTPs*iSpot;
                % get the desired measure of intensity
                switch opts.measure
                    case 'mean'
                        spotInts(rnge,iChan) = thisSpotInt(iSpot).intensity(:);
                    case 'max'
                        spotInts(rnge,iChan) = thisSpotInt(iSpot).intensity_max(:);
                    case 'median'
                        spotInts(rnge,iChan) = thisSpotInt(iSpot).intensity_median(:);
                    case 'min'
                        spotInts(rnge,iChan) = thisSpotInt(iSpot).intensity_min(:);
                end
                
            end
        end
        
        % correct spotInts for background
        bg = repmat(bg,nSpots,1);
        spotInts = spotInts - bg;
        
        % store cell's data
        allSpotInts = [allSpotInts; spotInts];
        allBg = [allBg; bg];
        % calculate cell-normalised data
        meanCtrl = nanmean(spotInts(:,opts.controlChan));
        allNorms = [allNorms; spotInts./meanCtrl];
        
    end
    
end

indNorms = allSpotInts./repmat(allSpotInts(:,opts.controlChan),1,3);    
    
%% Save and output data

% store data in analysis structure
analysis.channels = sort(opts.channels);
analysis.raw = allSpotInts;
analysis.bg = allBg;
analysis.norm.cellwise = allNorms;
analysis.norm.spotwise = indNorms;

% update stats - make a table of the statistics


% give results to output arguments
if nargout == 1
    varargout{1} = analysis;
end


end