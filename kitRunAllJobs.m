function kitRunAllJobs(jobsets)
% KITRUNALLJOBS Runs all jobs across multiple jobsets.
%
% Accepts cell array of filenames or cell array of jobset structs.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

% Ensure input is cell array.
if ~iscell(jobsets)
  jobsets = {jobsets};
end

nJobsets = length(jobsets);
kitLog('Running %d jobs',nJobsets);
for i = 1:nJobsets
  jobset = jobsets{i};
  if ischar(jobset)
    jobset = kitLoadJobset(jobset);
  end
  kitLog('Running job %i: %s',i,jobset.filename);

  % Run it.
  kitRunJob(jobset);
end
