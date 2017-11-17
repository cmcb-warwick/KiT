function jobset=kitConcatJobset(jobset1,jobset2)
% KITCONCATJOBSET Concatenates jobsets
%
%  Uses filename, options etc from first jobset.

% Sanity checks.
if ~strcmp(jobset1.movieDirectory,jobset2.movieDirectory)
  error('Different movie directorys');
end

% Concatenate movieFiles and crops.
jobset = jobset1;
concatFields = {'movieFiles','crop','cropSize'};
for i=1:length(concatFields)
  f = concatFields{i};
  if size(jobset.(f),2) == 1
    % In columns.
    jobset.(f) = [jobset.(f); jobset2.(f)];
  else
    % In rows.
    jobset.(f) = [jobset.(f) jobset2.(f)];
  end
end



