function job=kitTrackMovie(job,tasks)
% KITTRACKMOVIE Generates tracks for a single movie
%
%    JOB = KITTRACKMOVIE(JOB) Generates tracks for a single movie described by
%    JOB. Populates cell array field .dataStruct with results for each channel.
%
% Copyright (c) 2013 Jonathan W. Armond

tstart = tic;

if nargin<2
  tasks = 1:7;
end
% 1: finding spots
% 2: fitting plane
% 3: tracking spots
% 4: grouping sisters
% 5: extracting tracks
% 6: updating classes
% 7: aligning
% 8: intensity

if any(ismember(tasks,[1 2 8]))
  % Open movie and read metadata.
  if isfield(job,'metadata')
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie),'valid',job.metadata{job.index});
    job = kitSaveJob(job);
  else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie));
    job = kitSaveJob(job);
  end
end

opts = job.options;
nChannels = job.metadata.nChannels;

% Check which channels to analyze.
for c = 1:nChannels
  if strcmp(opts.coordMode{c}, 'none')
    ci(c) = 0;
  else
    ci(c) = 1;
  end
end
channels = find(ci);
if isempty(channels)
  error('No channels selected for analysis (movie has %d)',nChannels);
end
job.analyzedChannels = channels;

% Make dataStructs and take into account any changed options.
for c = channels
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
  % Update classes.
  for c = channels
    kitLog('Update kt classes in channel %d', c);
    job.dataStruct{c} = kitUpdateClass(job.dataStruct{c});
  end
  job = kitSaveJob(job);
end

if ismember(7,tasks)
  % Get alignment.
  for c = channels
    kitLog('Compute alignment in channel %d',c);
    job.dataStruct{c} = kitAlignFrames(job.dataStruct{c});
  end
  job = kitSaveJob(job);
end

if ismember(8,tasks)
  % Read spot intensity.
  for c = channels
    kitLog('Measure particle intensity in channel %d',c);
    job = kitLocalIntensityTracks(job, reader, job.metadata, c);
  end
  job = kitSaveJob(job);
end

% Gather diagnostics.
elapsed = toc(tstart);
for c = channels
  kitLog('Gather diagnostics in channel %d',c);
  job = kitDiagnostics(job,c,elapsed);
  fprintf('Diagnostics for channel %d\n', c);
  fprintf('-------------------------\n', c);
  kitPrintDiagnostics(job.dataStruct{c});
  fprintf('-------------------------\n\n', c);
end
job = kitSaveJob(job);

if any(ismember(tasks,[1 2 8]))
  reader.close();
  clear reader;
end

