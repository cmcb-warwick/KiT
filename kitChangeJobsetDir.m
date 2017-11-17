function jobset=kitChangeJobsetDir(jobset,newpath)
% KITCHANGEJOBSETDIR Change root directory of jobset
%
%    JOBSET = KITCHANGEJOBSETDIR(JOBSET,NEWPATH) Changes root directory of
%    jobset. The root directory is where the jobset mat-file is stored and
%    serves as the root directory for movie paths.
%
%    NEWPATH Optional. Specify the new path. If empty, dialog box is presented.
%
% Copyright (c) 2013 Jonathan W. Armond

oldpath = jobset.movieDirectory;

if nargin < 2
  newpath = uigetdir(oldpath,'Choose new root directory');
  if isemtpy(newpath) || newpath == 0
    error('Cancelled');
  end
end

jobset.movieDirectory = newpath;

% Replace in jobset filename.
jobset.filename = strrep(jobset.filename,oldpath,newpath);
