function job=kitTrackMovie(job,tasks)
% KITTRACKMOVIE Generates tracks for a single movie
%
%    JOB = KITTRACKMOVIE(JOB) Generates tracks for a single movie described by
%    JOB. Populates cell array field .dataStruct with results for each channel.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2016 C. A. Smith

tstart = tic;

if nargin<2
  switch job.options.jobProcess
    case 'zandt'
      tasks = 1:8;
    case 'zonly'
      tasks = [1,2,6];
    case 'chrshift'
      tasks = [1,6];
  end
  if ~any(job.options.intensity.execute)
    tasks = setdiff(tasks,9);
  end
end
% 1: finding spots
% 2: fitting plane
% 3: tracking spots
% 4: grouping sisters
% 5: extracting tracks
% 6: finding neighbour spots
% 7: updating classes
% 8: aligning
% 9: intensity

if any(ismember(tasks,[1 2 9]))
  % Open movie and read metadata.
  if isfield(job,'metadata')
    if iscell(job.metadata)
      job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);
  else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
  end
  job = kitSaveJob(job);
end

opts = job.options;
nChannels = min(job.metadata.nChannels,3); % enforce for now that there are
                                           % only 3 channels, as KiT can
                                           % only handle this many

% Check which channels to analyze.
for c = 1:nChannels
  if strcmp(opts.coordMode{c},'none')
    ci(c) = 0;
  elseif strcmp(opts.spotMode{c},'neighbour')
    ci(c) = -1;  
  else  
    ci(c) = 1;
  end
end
channels = find(ci==1);
neighChans = find(ci==-1);
if isempty(channels)
  error('No channels selected for analysis (movie has %d)',nChannels);
end
job.analyzedChannels = sort([channels neighChans]);

% Make dataStructs and take into account any changed options.
for c = [channels neighChans]
  ds = kitMakeMakiDatastruct(job,c);
  job.dataStruct{c}.dataProperties = ds.dataProperties;
end

if ismember(1,tasks)
  % Find 3D spot coordinates per frame.
  for c = channels
    kitLog('Finding particle coordinates in channel %d',c);
    job = kitFindCoords(job, reader, c);
    job = kitSaveJob(job);
    % Give up if spot finding failed.
    if isfield(job.dataStruct{c},'failed') && job.dataStruct{c}.failed
      warning('Giving up on job.');
      return
    end
  end
end

if ismember(2,tasks)
  % Fit plane in chosen channel.
  planeChan = job.options.coordSystemChannel;
  kitLog('Fitting plane in channel %d', planeChan);
  if strcmp(opts.coordMode{planeChan}, 'none')
    % No spot tracking in plane channel so populate dataStruct.
    job.dataStruct{planeChan} = kitMakeMakiDatastruct(job, planeChan);
  end
  if strcmp(job.options.coordSystem, 'register')
    job.dataStruct{planeChan} = kitRegisterFrames(job,reader,job.dataStruct{planeChan});
  else
    job.dataStruct{planeChan} = kitFitPlane(job,reader,job.dataStruct{planeChan},planeChan,0);
  end
  % Transform other channels.
  for c = channels
    if c ~= planeChan
      kitLog('Transforming coordinates in channel %d to plane from channel %d',c,planeChan);
      job.dataStruct{c}.planeFit = job.dataStruct{planeChan}.planeFit;
      job.dataStruct{c} = kitFitPlane(job,reader,job.dataStruct{c},c,1);
    end
  end
  job = kitSaveJob(job);
else
  for c = channels
    if ~isfield(job.dataStruct{c},'planeFit')
      job.dataStruct{c}.planeFit = [];
    end
  end
end

if ismember(3,tasks)
  % Track spots.
  for c = channels
    kitLog('Tracking particles in channel %d', c);
    job.dataStruct{c} = kitGenerateTracks(job.dataStruct{c});
  end
  job = kitSaveJob(job);
end

if ismember(4,tasks)
  % Group sisters.
  for c = channels
    kitLog('Grouping sisters in channel %d', c);
    job.dataStruct{c} = kitGroupSisters(job.dataStruct{c},opts.debug.groupSisters);
  end
  job = kitSaveJob(job);
end

if ismember(5,tasks)
  % Extract individual tracks.
  for c = channels
    if isfield(job.dataStruct{c},'planeFit')
      kitLog('Extracting individual tracks in channel %d', c);
      job = kitExtractTracks(job, c);
    end
  end
  job = kitSaveJob(job);
end

if ismember(6,tasks)
  % Find neighbouring 3D spot coordinates per frame.
  for c = neighChans
    if ismember(1,tasks)
      kitLog('Finding particle coordinates in channel %d',c);
      job = kitFindCoords(job, reader, c);
    end
    if ismember(2,tasks)
      % Transform coordinates into plane.
      kitLog('Transforming coordinates in channel %d to plane from channel %d',c,planeChan);
      planeChan = job.options.coordSystemChannel;
      job.dataStruct{c}.planeFit = job.dataStruct{planeChan}.planeFit;
      job.dataStruct{c} = kitFitPlane(job,reader,job.dataStruct{c},c,1);
    end
    if ismember(3,tasks)
        % Assemble tracks and sisterList from initCoord and planeFit structures
        % in neighbour channel.
        if strcmp(opts.jobProcess,{'zandt'})
          job.dataStruct{c} = kitAssembleNeighbourStructs(job.dataStruct,c,opts);
          % Extract tracks from tracks and sisterList.
          kitLog('Extracting individual tracks in channel %d', c);
          job = kitExtractTracks(job, c);
        end
    end
  end
  job = kitSaveJob(job);
end

if ismember(7,tasks)
  % Update classes.
  for c = [channels neighChans]
    kitLog('Update kt classes in channel %d', c);
    job.dataStruct{c} = kitUpdateClass(job.dataStruct{c});
  end
  job = kitSaveJob(job);
end

if ismember(8,tasks)
  % Get alignment.
  for c = [channels neighChans]
    kitLog('Compute alignment in channel %d',c);
    job.dataStruct{c} = kitAlignFrames(job.dataStruct{c});
  end
  job = kitSaveJob(job);
end

if ismember(9,tasks)
  % Read spot intensity.
  intChans = find(job.options.intensity.execute);
  for c = intChans
    kitLog('Measure particle intensity in channel %d',c);
    job = kitLocalIntensity(job, reader, job.metadata, c, opts.intensity);
  end
  job = kitSaveJob(job);
end

% Gather diagnostics.
elapsed = toc(tstart);
for c = sort([channels neighChans])
  kitLog('Gather diagnostics in channel %d',c);
  job = kitDiagnostics(job,c,elapsed);
  fprintf('Diagnostics for channel %d\n', c);
  fprintf('-------------------------\n', c);
  kitPrintDiagnostics(job.dataStruct{c},opts.jobProcess);
  fprintf('-------------------------\n', c);
end
fprintf('\n')
job = kitSaveJob(job);

if any(ismember(tasks,[1 2 9]))
  reader.close();
  clear reader;
end

