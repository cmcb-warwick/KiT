function jobs = chrsFilterSpots(jobs,varargin)
% CHRSFILTERSPOTS Filters spots located for a chromatic shift job to allow
%    more accurate calculation of chromatic shift.
%
%    jobs = CHRFILTERSPOTS(JOBS,...) Returns a JOB structure, typically for
%         chromatic shift calculation, with initCoord coordinates filtered
%         for location in the imaging plane, distance from its nearest
%         neighbour, spot brightness for a given channel, or a combination
%         of the three.
%
%    Options, defaults in {}:-
%
%    appliedChans: [1 2] or subset of [1:3]. Channels for which filters
%         apply.
%
%    intensityFilter: {0} or number between 0 and 1. The minimum ratio of a
%         given spot's intensity divided by the average of the brightest
%         5% of all spots.
%         NB. If less that 20 spots are localised, the ratio is calculated
%         relative to the brightest spot.
%
%    neighbourFilter: {0} or distance in microns. The minimum distance
%         allowed between two spots.
%
%    referenceChan: {1} or number between 1 and 3. The channel on which
%         filtering is based.
%
%    regionFilter: {1} or an odd number. The number, n, of divisions over
%         the xy-plane over which to filter data to only the central n'th
%         in each the x- and y-coordinates.
%         NB. A value of 1 results in no filter applied.
%         NB2. If an even number is provided, it will be decreased by 1.
%
%    revert: {1} or 0. Whether or not to revert to the raw initCoord
%         structure in the event that filtering has been performed
%         previously.
%
%
% Copyright (C) 2015 C. A. Smith

if nargin<1 || isempty(jobs)
    error('Please provide a chromatic shift job.')
end

% set default options
opts.appliedChans = [1 2];
opts.intensityFilter = 0;
opts.neighbourFilter = 0;
opts.referenceChan = 1;
opts.regionFilter = 1;
opts.revert = 1;
% process user-defined options
opts = processOptions(opts,varargin{:});

for iJob = 1:length(jobs)

    job = jobs{iJob};
    
    % reset some structures
    coords = [];
    amps = [];
    
%% Data handling

    % get dataStructs from reference channel
    refDataStruct = job.dataStruct{opts.referenceChan};

    % get original initCoords for all channels being analysed, and produce new
    % initCoord structures
    for iChan = unique([opts.referenceChan opts.appliedChans])

      % get the raw initCoord structure if present, and revert is active
      if isfield(job.dataStruct{iChan},'failed') && ~job.dataStruct{iChan}.failed
        if opts.revert && isfield(job.dataStruct{iChan},'rawInitCoord')
          origInitCoord{iChan} = job.dataStruct{iChan}.rawInitCoord;
        else
          origInitCoord{iChan} = job.dataStruct{iChan}.initCoord(1);
        end
      else
        return
      end
      % copy this information to a new structure for later
      newInitCoord{iChan} = origInitCoord{iChan};

    end

    % get the number of spots, and check that this is equal in all channels
    nSpots = size(origInitCoord{opts.referenceChan}.allCoord,1);
    for iChan = setdiff(opts.appliedChans,opts.referenceChan)
      tempNspots = size(origInitCoord{iChan}.allCoord,1);
      if tempNspots ~= nSpots
        error('Have a different number of points in the channel %i compared to channel %i.',iChan,opts.referenceChan)
      end
    end

    % get pixelSize and imageSize information
    pixelSize = job.metadata.pixelSize;
    if isempty(job.ROI)
      imageSize = job.metadata.frameSize;
    else
      imageSize = job.ROI.cropSize;
    end

    % get coordinate and intensity data
    for iChan = unique([opts.referenceChan opts.appliedChans])
      coords(:,:,iChan) = origInitCoord{iChan}.allCoord(:,1:3);
      amps(:,iChan) = origInitCoord{iChan}.amp(:,1);
    end

%% Find which data requires filtering

    % region filtering
    nRegions = opts.regionFilter;
    % designate a vector of ones to accept all current spots
    regionPass = ones(nSpots,1);
    filterParams.numberRegions = [];

    if nRegions < 1

      % ensure that a sensible number of regions is given
      warning('Number of regions given for filtering must be larger than or equal to 1. Region filtering aborted.')

    elseif all(nRegions~=[1 2])

      % check that the number of regions is odd
      if mod(nRegions,2) == 0
        warning('Number of regions given for filtering must be odd. Number of regions decreased from %i to %i.',nRegions,nRegions-1)
        nRegions = nRegions - 1;
      end

      % find the location of the central division lines, round to the nearest
      % pixel, then convert to microns
      regionDivisions(1,:) = (nRegions-1)/(2*nRegions)*imageSize(1:2);
      regionDivisions(2,:) = (nRegions+1)/(2*nRegions)*imageSize(1:2);
      regionDivisions = round(regionDivisions);
      regionDivisions = regionDivisions.*repmat(pixelSize(1:2),2,1);

      % find which spots are within the central region, saving the filter
      % parameters
      for iSpot = 1:nSpots
          regionPass(iSpot) = all(coords(iSpot,1:2,opts.referenceChan)>regionDivisions(1,:));
          regionPass(iSpot) = all(coords(iSpot,1:2,opts.referenceChan)<regionDivisions(2,:));
      end
      filterParams.numberRegions = nRegions;

    end

    % intensity filtering
    minIntensity = opts.intensityFilter;
    % designate a vector of ones to accept all current spots
    intensityPass = ones(nSpots,1);
    filterParams.minIntensity = [];

    if minIntensity > 1

      % ensure that minimum intensity ratio given is smaller than or equal to 1  
      warning('Minimum intensity ratio cannot be larger than 1. Intensity filtering aborted.')

    else

      % sort the intensities in decending order
      maxAmps = sort(amps(:,opts.referenceChan),1,'descend');
      % get 20% of total number of spots, and find the average of this many
      % spots
      nAverage = ceil(nSpots/5);
      averageMaxInt = nanmean(maxAmps(1:nAverage));
      % find the minumum intensity required compared to this maximum
      % measurement
      minIntensity = minIntensity*averageMaxInt;

      % find which spots are within the intensity limitations, saving the
      % filter parameters
      intensityPass = amps(:,opts.referenceChan)>minIntensity;
      filterParams.minIntensity = minIntensity;

    end

    % neighbour filter
    minNNdist = opts.neighbourFilter;
    % designate a vector of ones to accept all current spots
    neighbourPass = ones(nSpots,1);
    filterParams.minNeighbourSeparation = [];

    if minNNdist > max(imageSize.*pixelSize)

      % check that minimum neighbour separations are within the image size
      warning('Minimum neighbour separation cannot be larger than the image size. Neighbour filtering aborted.')

    elseif minNNdist > 0

      % get the distances between spots and its nearest neighbour
      nnDists = createDistanceMatrix(coords(:,:,opts.referenceChan),coords(:,:,opts.referenceChan));
      nnDists = sort(nnDists,2);
      nnDists = nnDists(:,2);

      % find which spots are outside the accepted neighbour separation, saving
      % the filter parameters
      neighbourPass = nnDists>minNNdist;
      filterParams.minNeighbourSeparation = minNNdist;

    end

    % find all coordinates that fail any criteria
    failed = (intensityPass+neighbourPass+regionPass) < 3;

    % check if any filtering has occurred, and if not then avoid production of
    % rawInitCoord
    if sum(failed) == 0
      fprintf('No filtering occurred under the parameters given. No new initCoord will be produced.\n')
      return
    end

    for iChan = opts.appliedChans

      % remove unwanted data from the new initCoord structures
      newInitCoord{iChan}.allCoord(failed,:) = NaN;
      newInitCoord{iChan}.allCoordPix(failed,:) = NaN;
      newInitCoord{iChan}.amp(failed,:) = NaN;
      newInitCoord{iChan}.bg(failed,:) = NaN;
      % include new nSpots value, and filter parameters
      newInitCoord{iChan}.nSpots = sum(1-failed);
      newInitCoord{iChan}.filterParams = filterParams;

      % back up the original initCoord structures
      if ~isfield(job.dataStruct{iChan},'rawInitCoord')
        job.dataStruct{iChan}.rawInitCoord = origInitCoord{iChan};
      end
      % store new initCoord information
      job.dataStruct{iChan}.initCoord = newInitCoord{iChan};

    end

    % back up and save job
    job = kitSaveJob(job);
    jobs{iJob} = job;

end

end