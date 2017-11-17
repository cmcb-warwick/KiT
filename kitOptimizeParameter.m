function diagnostics=kitOptimizeParameter(jobset,jobId,parameter,values,varargin)
% KITOPTIMIZEPARAMETER Run tracking over sweep of parameter for optimization
%
%    DIAGNOSTICS = KITOPTIMIZEPARAMETER(JOBSET,JOBID,PARAMETER,VALUES,...)
%    Run tracking for job numbered JOBID from JOBSET, over sweep of PARAMETER
%    across VALUES for optimization. PARAMETER should be a string representing a
%    valid field name of JOB.OPTIONS. VALUES should be a vector containing the
%    values for PARAMETER, or a matrix where each row is a value for PARAMETER.
%
%    Optional parameters, string/value pairs.
%
%    'channel':  Specify tracking channel, default = 1.
%    'parallel': Run jobs in batch mode, default = 0.
%    'tasks': Tracking phases to run, default = 1:8.
%    'existing': Use existing tracking results by loading job first, default = 0.
%
% Copyright (c) 2013 Jonathan W. Armond

opts.channel = 1;
opts.parallel = 0;
opts.tasks = 1:8;
opts.existing = 0;
opts = processOptions(opts,varargin{:});

if ~ischar(parameter) || ~isfield(jobset.options,parameter)
  error('PARAMETER should be a string representing field of JOB.OPTIONS');
end

% Ensure column vector, if vector.
if isvector(values) && size(values,2) > 1
  values = values';
end

% Upgrade jobset options, if required.
if ~isfield(jobset.options,'jobsetVersion') || ...
    jobset.options.jobsetVersion < kitVersion(2)
  defJob = kitDefaultOptions();
  jobset = structCopyMissingFields(jobset,defJob);
end

nValues = size(values,1);
for i = 1:nValues
  kitLog('Setting option %s to %s',parameter,num2str(values(i),'%g '));

  if opts.existing
    job = kitLoadJob(jobset,jobId);
  else
    % Create job.
    job = jobset;
    job.movie = jobset.movieFiles{jobId};
    job.index = jobId;
    job.crop = jobset.crop{jobId};
    job.cropSize = jobset.cropSize{jobId};
    % Update versions, may be different to jobset creator.
    job.version = kitVersion();
    job.matlabVersion = version;
  end
  job.options.(parameter) = values(i,:);
  % Inhibit saving, since all jobs would use same file.
  job.options.disableSave = 1;

  % Submit tracking job.
  if opts.parallel
    kitLog('Submitting tracking job %d', i);
    batchJob(i) = batch(@kitTrackMovie, 1, {job, opts.tasks});
  else
    kitLog('Running tracking job %d', i);
    job = kitTrackMovie(job, opts.tasks);
    diagnostics(i) = job.dataStruct{opts.channel}.diagnostics;
  end
end

% Wait for tracking jobs.
if opts.parallel
  for i = 1:nValues
    kitLog('Waiting for %d tracking jobs to complete', nValues-i+1);
    b = batchJob(i);
    wait(b);
    diary(b);

    if strcmp(b.State,'finished')
      % Get result.
      r = fetchOutputs(b);
      job = r{1};
      diagnostics(i) = job.dataStruct{opts.channel}.diagnostics;
    else
      kitLog('Job %d failed', i);
    end

    delete(b);
  end
end

% Save output.
[~,jobsetName] = fileparts(job.filename);
[moviePath,movieName] = fileparts(job.movie);
outputName = fullfile(job.movieDirectory,moviePath,...
                      ['diagnostics-' jobsetName '-' parameter '.mat']);
save(outputName,'diagnostics','parameter','values');
