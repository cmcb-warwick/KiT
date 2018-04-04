function [chrShift,nDataPoints] = chrsCalculateChromaticShift(data,chanVect,varargin)
% CHRSCALCULATECHROMATICSHIFT Produces a vector for chromatic shift from a
%    set of images previously processed for chromatic shift measurement.
%
%    [chrShift,nDataPoints] = KITCALCULATECHROMATICSHIFT(JOBS,CHANVECT,...)
%    Chromatic shift is calculated from a structure of JOBS with spots
%    located in each of the two channels given in the channel vector,
%    CHANVECT, which also defines the direction in which chromatic shift is
%    calculated.
%
%    Options, defaults in {}:-
%
%    filtered: 0 or {1}. Whether or not to use filtered data.
%
%    interphaseToMetaphase: {0} or 1. Whether or not to impose a z-
%         directional adjustment (+78.2nm) when using interphase cells to
%         calculate chromatic shift for use in metaphase cells.
%         N.B. This is only required when using the dual camera system on
%         the McAinsh ultraView spinning disk microscope.
%
%
% Copyright (C) 2017 C. A. Smith

% set default options
opts.filtered = 1;
opts.interphaseToMetaphase = 0;
% process options
opts = processOptions(opts,varargin{:});

% if data is a zetaStructure, provide results from this
if isstruct(data) && isfield(data,'chrsVersion')
    
    kitLog('Pre-processed zetaStructure detected.')
    
    chrShift = zeros(1,6);
    if opts.filtered
        % find median and standard deviation of chromatic shift measurements
        chrShift(1) = nanmedian(data.filtered.zeta.x(:));
        chrShift(2) = nanmedian(data.filtered.zeta.y(:));
        chrShift(3) = nanmedian(data.filtered.zeta.z(:));
        chrShift(4) = nanstd(data.filtered.zeta.x(:));
        chrShift(5) = nanstd(data.filtered.zeta.y(:));
        chrShift(6) = nanstd(data.filtered.zeta.z(:));
        nDataPoints = sum(~isnan(data.filtered.zeta.x));
    else
        % find median and standard deviation of chromatic shift measurements
        chrShift(1) = nanmedian(data.raw.zeta.x(:));
        chrShift(2) = nanmedian(data.raw.zeta.y(:));
        chrShift(3) = nanmedian(data.raw.zeta.z(:));
        chrShift(4) = nanstd(data.raw.zeta.x(:));
        chrShift(5) = nanstd(data.raw.zeta.y(:));
        chrShift(6) = nanstd(data.raw.zeta.z(:));
        nDataPoints = sum(~isnan(data.raw.zeta.x));
    end
    % give a warning in the case that there are less than 20 datapoints
    if nDataPoints < 20
        warning('Have very few datapoints across all chromatic shift movies. Consider revising.')
    end
    %correct z-coordinate from interphase to metaphase if necessary
    if opts.interphaseToMetaphase
        chrShift(3) = chrShift(3) + 0.0782; % correct for McA spinning disk dual camera
    end
    
    kitLog('Chromatic shift calculated: [%.1f %.1f %.1f] ± [%.1f %.1f %.1f] (n = %i)',...
        chrShift(1)*1000,chrShift(2)*1000,chrShift(3)*1000,chrShift(4)*1000,chrShift(5)*1000,chrShift(6)*1000,nDataPoints)
    return
    
end

% check that a channel vector has been provided
if nargin<2 || isempty(chanVect)
    chanVect = [1 2];
end

% find number of jobs
nJobs = length(data);
% preallocate vector for accumulating chromatic shift measurements, zeta
accumZeta = [];

% loop over jobs
for iJob = 1:nJobs
    % get number of channels
    nChans = length(data{iJob}.dataStruct);
    % check that this is enough based on the provided channel vector
    if any(chanVect>nChans)
        chrShift = zeros(1,6);
        nDataPoints = NaN;
        return
    end
    
    % get coordinates, calculate zeta
    if isfield(data{iJob}.dataStruct{chanVect(2)},'initCoord') && ~data{iJob}.dataStruct{chanVect(2)}.failed
      
      % get coordinates for each channel
      if opts.filtered
        coords1 = data{iJob}.dataStruct{chanVect(1)}.initCoord.allCoord(:,1:3);
        coords2 = data{iJob}.dataStruct{chanVect(2)}.initCoord.allCoord(:,1:3);
      else
        % check for a rawInitCoord in each channel dataStruct
        % channel 1
        if isfield(data{iJob}.dataStruct{chanVect(1)},'rawInitCoord')
            coords1 = data{iJob}.dataStruct{chanVect(1)}.rawInitCoord.allCoord(:,1:3);
        else
            coords1 = data{iJob}.dataStruct{chanVect(1)}.initCoord.allCoord(:,1:3);
        end
        % channel 2
        if isfield(data{iJob}.dataStruct{chanVect(2)},'rawInitCoord')
            coords2 = data{iJob}.dataStruct{chanVect(2)}.rawInitCoord.allCoord(:,1:3);
        else
            coords2 = data{iJob}.dataStruct{chanVect(2)}.initCoord.allCoord(:,1:3);
        end
      end
      
      % calculate distances and accumulate for all experiments
      zeta = coords2 - coords1;
      accumZeta = [accumZeta; zeta];
    end
    
end

% find number of datapoints accumulated across the entire jobset
nDataPoints = sum(~isnan(prod(accumZeta,2)));
if nDataPoints == 0
    chrShift = zeros(1,6);
    warning('No data. Setting chromatic shift = [0 0 0].')
else
    % give a warning in the case that there are less than 20 datapoints
    if nDataPoints < 20
        warning('Have very few datapoints across all chromatic shift movies. Consider revising.')
    end
    % find median and standard deviation of chromatic shift measurements
    chrShift = nanmedian(accumZeta);
    chrShift(4:6) = nanstd(accumZeta);
    %correct z-coordinate from interphase to metaphase if necessary
    if opts.interphaseToMetaphase
        chrShift(3) = chrShift(3) + 0.0782; % correct for McA spinning disk dual camera
    end
end

kitLog('Chromatic shift calculated: [%.1f %.1f %.1f] ± [%.1f %.1f %.1f] (n = %i)',...
        chrShift(1)*1000,chrShift(2)*1000,chrShift(3)*1000,chrShift(4)*1000,chrShift(5)*1000,chrShift(6)*1000,nDataPoints)

end