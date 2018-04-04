function jobs=kitLoadAllJobs(jobset,subset)
% KITLOADALLJOBS Loads all jobs specified in a jobset into cell array.
%
%    KITLOADALLJOBS(JOBSET,SUBSET) Loads a SUBSET of jobs specified in the
%    JOBSET into a cell array. If a SUBSET is not provided, all jobs within
%    the JOBSET are loaded.
%
% Created by: J. W. Armond
% Copyright (c) 2016 C. A. Smith
    
if isfield(jobset,'exclude')
  exclude = jobset.exclude;
else
  exclude = [];
end

n = length(jobset.ROI);
if nargin<2 || isempty(subset)
  subset = 1:n;
elseif ~isempty(setdiff(subset,1:n))
  error('Subset too large for the jobset provided.')
end

jobs = {};
for i=subset
  % Exclude cells.
  if ismember(i,exclude)
    continue
  end

  jobs{end+1} = kitLoadJob(jobset,i);
end
