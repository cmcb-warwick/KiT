function job = kitSaveJob(job)
% KITSAVEJOB Saves job struct with standardized name
%
%    JOB = KITSAVEJOB(JOB) Saves job struct as mat-file with name based on movie
%    filename: kittracking_jobsetFilename_movieFilename.mat
%
% Copyright (c) 2013 Jonathan W. Armond

% Certain modes prohibit saving of output.
if (isfield(job.options,'disableSave') && job.options.disableSave) || ...
    (isfield(job,'disableSave') && job.disableSave)
  return;
end

% Generate output name.
outputName = kitGenerateOutputFilename(job);

if ~isfield(job,'output')
    % Generate output name.
    job.output = outputName;
end

% Save mat.
save(job.output, '-struct', 'job', '-v7.3');
