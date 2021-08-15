function job = kitSaveJob(job)
% KITSAVEJOB Saves job struct with standardized name
%
%    JOB = KITSAVEJOB(JOB) Saves job struct as mat-file with name based on movie
%    filename: kittracking_jobsetFilename_movieFilename.mat
%
% Copyright (c) 2019 Jonathan U. Harrison and Jonathan W. Armond

% Certain modes prohibit saving of output.
if (isfield(job.options.debug,'disableSave') && job.options.debug.disableSave) || ...
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
try
  save(job.output, '-struct', 'job', '-v7.3');
catch 
  %had some errors using shared server relating to permissions
  %instead try to save in temporary location and copy across to movie dir
  warning('Initially uanble to save file; saving in temporary location and copying across');
  savename = tempname(); %get a temporary location to put file instead
  save(savename,'-struct','job','-v7.3');
  system(sprintf('cp %s.mat %s',savename,strrep(job.output,' ', '\ ')),'-echo'); %copy across, be careful about spaces
  system(sprintf('rm %s.mat',savename),'-echo'); %clean up afterwards
end

