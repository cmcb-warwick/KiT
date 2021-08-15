function kitSaveJobset(jobset)
% KITSAVEJOBSET Save jobset in mat-file
%
% Copyright (c) 2019 Jonathan U. Harrison and Jonathan W. Armond
try
  save(jobset.filename,'-struct','jobset','-v7.3');
catch
  %had some errors using shared server relating to permissions
  %instead try to save in temporary location and copy across to movie dir
  warning('Initially uanble to save file; saving in a temporary location and copying across');
  savename = tempname(); %get a temporary location to put file instead
  save(savename,'-struct','jobset','-v7.3');
  system(sprintf('cp %s.mat %s',savename,strrep(jobset.filename,' ', '\ ')),'-echo'); %copy across, be careful about spaces in filenames
  system(sprintf('rm %s.mat',savename),'-echo'); %clean up afterwards
end

