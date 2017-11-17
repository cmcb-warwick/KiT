function analysis=cuplLoadData(analysis,jobset,channel)
% CUPLLOADDATA Load dataset into an struct
%
% Modified to integrate into KiT.
%
% Copyright (c) 2013 Jonathan Armond

% Allocate and initialise storage.
an = analysis;
an.version = kitVersion();
nROIs = length(jobset.ROI);
sisterCount = zeros(nROIs,1); % Nx1 (N: # of cells)
trackCount = zeros(nROIs,1); % Nx1 (N: # of cells)
an.sisterCellIdx = []; % 1xN (N: # of sisters)
an.sisterCoords1 = []; % Tx3N (T: time points, N: # of sisters)
an.sisterCoords2 = [];
an.trackCellIdx = []; % 1xN (N: # of tracks)
an.trackCoords = []; % Tx3N (T: time points, N: # of tracks)
an.trackInt = [];
an.hasTracks = 0;
an.hasTrackInt = 0;

% Load cell data.
for i=1:nROIs
  try
    job = kitLoadJob(jobset, i);
  catch me
    if strcmp(me.identifier,'MATLAB:load:couldNotReadFile')
      warning('File missing for job %d.',i)
    end
  end
  empty = 0;
  if ~isfield(job,'dataStruct') || length(job.dataStruct) < channel
    empty = 1;
  else
    data = job.dataStruct{channel};
    if isempty(data)
      empty = 1;
    end
  end
  if empty
    warning('Empty data in job %d. Did you choose the right channel?',i);
    continue;
  end

  if isfield(data,'sisterList') && ~isempty(data.sisterList) && ~isempty(data.sisterList(1).coords1)
    % Sister track coords, stored as sequential coordinates, e.g.
    % an.sisterCoords1(:,1:3) selects x,y,z of first sister.
    nSisters = length(data.sisterList);
    sisterCount(i) = nSisters;
    an.sisterCellIdx = [an.sisterCellIdx repmat(i,1,nSisters)];

    % Expand if necessary.
    nFrames = size(data.sisterList(1).coords1,1);
    an.sisterCoords1 = expand(an.sisterCoords1, nFrames);
    an.sisterCoords2 = expand(an.sisterCoords2, nFrames);

    % Concatenate data.
    for j=1:nSisters
      an.sisterCoords1(1:nFrames,end+1:end+3) = data.sisterList(j).coords1(:,1:3);
      an.sisterCoords2(1:nFrames,end+1:end+3) = data.sisterList(j).coords2(:,1:3);
    end
  end

  if isfield(data,'trackList') && ~isempty(data.trackList)
    % Track coords. Not the same as sister track coords as it includes parts
    % of tracks unable to be assigned a pair.
    nTracks = length(data.trackList);
    trackCount(i) = nTracks;
    an.trackCellIdx = [an.trackCellIdx repmat(i,1,nTracks)];

    % Expand if necessary.
    nFrames = size(data.trackList(1).coords,1);
    an.trackCoords = expand(an.trackCoords, nFrames);

    % Concatenate data.
    for j=1:nTracks
      an.trackCoords(1:nFrames,end+1:end+3) = data.trackList(j).coords(:,1:3);
    end
    an.hasTracks = 1;
  end

  if isfield(data,'trackInt') && ~isempty(data.trackInt)
    % Track intensities.
    nTracks = length(data.trackInt);
    if ~isfield(data,'trackList') || nTracks ~= length(data.trackList)
      warning('Size of trackInt does not match size of trackList');
    end

    if isempty(an.trackInt)
      % Create struct array for number of channels.
      nChannels = size(data.trackInt(1).intensity,2);
      trackInt(1:nChannels) = struct('mean',[],'max',[],'min',[],'median',[]);
      an.trackInt = trackInt; % KLUDGE can't reassign to struct
    else
      if nChannels ~= length(an.trackInt)
        error(['Number of intensity channels inconsistent. Now I don''t know ' ...
               'which channel goes where...']);
      end
    end

    for j=1:nChannels
      % Expand if necessary.
      nFrames = size(data.trackInt(1).intensity,1);
      an.trackInt(j).mean = expand(an.trackInt(j).mean, nFrames);
      an.trackInt(j).max = expand(an.trackInt(j).max, nFrames);
      an.trackInt(j).min = expand(an.trackInt(j).min, nFrames);
      an.trackInt(j).median = expand(an.trackInt(j).median, nFrames);
    end

    for j=1:nTracks
      for k=1:nChannels
        nFrames = size(data.trackInt(1).intensity,1);
        an.trackInt(k).mean(1:nFrames,end+1) = data.trackInt(j).intensity(:,k);
        an.trackInt(k).max(1:nFrames,end+1) = data.trackInt(j).intensity_max(:,k);
        an.trackInt(k).min(1:nFrames,end+1) = data.trackInt(j).intensity_min(:,k);
        an.trackInt(k).median(1:nFrames,end+1) = data.trackInt(j).intensity_median(:,k);
      end
    end
    an.hasTrackInt = 1;
  end
end

an.loadedSisterCount = sisterCount;
an.loadedTrackCount = trackCount;
an.nLoadedSisters = sum(sisterCount);
an.nLoadedTracks = sum(trackCount);
an.nLoadedCells = nROIs;
an.maxTrackLength = size(an.sisterCoords1,1);
% Time vector.
t = job.metadata.frameTime;
an.dt = t(1,2)-t(1,1);
an.time = (0:an.maxTrackLength-1)'*an.dt;
an.nChannels = length(an.trackInt);

% Record stage.
an.stages = union(an.stages,'loaded');

% Unalias analysis.
analysis = an;

%% LOCAL FUNCTIONS

function m=expand(m,rows)
% Expand m to have at least rows.
  sizeDiff = nFrames - size(m,1);
if ~isempty(m) && sizeDiff > 0
  m = [m; nan(sizeDiff,size(m,2))];
end
end

end
