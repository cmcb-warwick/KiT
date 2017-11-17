function jobs=kitLoadAllJobs(jobset)
% KITLOADALLJOBS Loads all jobs specified in a jobset into cell array.
%
% Copyright (c) 2013 Jonathan Armond

if isfield(jobset,'exclude')
  exclude = jobset.exclude;
else
  exclude = [];
end

n = length(jobset.ROI);
jobs = {};
for i=1:n
  % Exclude cells.
  if ismember(i,exclude)
    continue
  end

  jobs{end+1} = kitLoadJob(jobset,i);
end
