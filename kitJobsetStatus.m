function kitJobsetStatus(jobset)
% KITJOBSETSTATUS Display completion stage of running jobs

nJobs = length(jobset.movieFiles);
status = cell(nJobs);
for i=1:nJobs
  job = kitLoadJob(jobset,i);

  % Progressively update status until field not found.
  if isfield(job,'metadata')
    status{i} = 'spot finding';
  else
    status{i} = 'not started or failed';
  end

  if ~isfield(job,'dataStruct')
    continue
  end

  % Find first non-empty dataStruct.
  l = cellfun(@(x) length(x), job.dataStruct);
  idx = find(l,1);
  if isempty(idx)
    continue
  end

  ds = job.dataStruct{idx};
  if isfield(ds,'initCoord')
    status{i} = 'plane fitting';
  else
    continue
  end

  if isfield(ds,'planeFit')
    status{i} = 'making tracks';
  else
    continue
  end

  if isfield(ds,'tracks')
    status{i} = 'grouping sisters';
  else
    continue
  end

  if isfield(ds,'sisterList')
    status{i} = 'extracting tracks';
  else
    continue
  end

  if isfield(ds,'trackList')
    status{i} = 'updating classes';
  else
    continue
  end

  if isfield(ds,'updatedClass')
    status{i} = 'measuring track intensity';
  else
    continue
  end

  if isfield(ds,'trackInt')
    status{i} = 'gathering diagnostics';
  else
    continue
  end

  if isfield(ds,'diagnostics')
    status{i} = 'finished';
  else
    continue
  end
end

for i=1:nJobs
  fprintf('Job %03d: %s\n',i,status{i});
end
