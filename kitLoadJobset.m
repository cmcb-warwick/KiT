function jobset=kitLoadJobset(filename)
% KITLOADJOBSET Load jobset into struct
%
%    JOBSET = KITLOADJOBSET(FILENAME) Loads jobset from mat-file FILENAME into
%    JOBSET struct. Root directory of JOBSET is automatically set to the
%    directory part of FILENAME. This makes it easy to move the jobset and
%    movies around together, and analyze them in different locations.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<1
  [filename,pathname] = uigetfile(...
    '*.mat','Locate jobset to load');

  if isequal(filename,0) || isequal(pathname,0)
    error('No jobset filename specified.');
  end

  filename = fullfile(pathname,filename);
end

if exist(filename,'file') == 0
  error('File not found: %s',filename);
end

% Get directory part of filename.
pathstr = fileparts(filename);

jobset = load(filename);

% Change root directory.
jobset = kitChangeJobsetDir(jobset, pathstr);
