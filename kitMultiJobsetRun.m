function kitMultiJobsetRun(jobsets,parallel)
% KITMULTIJOBSETRUN Runs all jobs from multiple jobsets in parallel or serially.
%
% Accepts cell array of filenames or cell array of jobset structs.
% Copyright (c) 2013 Jonathan W. Armond

if nargin<2
  parallel = 0;
end

% Ensure input is cell array.
if ~iscell(jobsets)
  jobsets = {jobsets};
end

nJobsets = length(jobsets);
kitLog('Running %d jobsets',nJobsets);
for i = 1:nJobsets
  jobset = jobsets{i};
  if ischar(jobset)
    jobset = kitLoadJobset(jobset);
  end
  kitLog('Running jobset: %s',jobset.filename);

  % Run it.
  kitRunJobs(jobset,'parallel',parallel);
end
