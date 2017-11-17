function outputName=kitGenerateOutputFilename(job)
% KITGENERATEOUTPUTFILENAME Generate filename for output mat based on jobset name

[~,jobsetName] = fileparts(job.filename);
[moviePath,movieName] = fileparts(job.movie);


if isfield(job,'jobsetVersion') && job.jobsetVersion >= 5
  % Add index near front of filename to improve filesystem sorting.
  fileName = ['kittracking' num2str(job.index,'%03d') '-' jobsetName '-' movieName];
else
  fileName = ['kittracking-' jobsetName '-' movieName];
end
if isfield(job,'variantName')
  % Add variant to name for testing purposes.
  fileName = [fileName '-' job.variantName];
end
if isfield(job,'jobsetVersion') && job.jobsetVersion > 2 && job.jobsetVersion < 5
  % Add index to filename.
  fileName = [fileName '-' num2str(job.index)];
end
fileName = [fileName '.mat'];

% Generate output name including jobset name.
outputName = fullfile(job.movieDirectory, moviePath, fileName);
