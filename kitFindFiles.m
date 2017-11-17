function filenames=kitFindFiles(directory,regexpr,incHidden,showprog)
% KITFINDFILESWITHEXTENSION Recursively search for files by regexp.
%
%   FILENAMES = KITFINDFILESWITHEXTENSION(DIRECTORY,REGEXPR,INCHIDDEN)
%
%   INPUT DIRECTORY: Directory to start recursive search.
%
%         REGEXPR: Regular expression for matching file names.
%
%         INCHIDDEN: Ignore files prefixed with '.' (0, default)
%                    or include them (1).
%
%         SHOWPROG: Shows progress bar when taking a long time.
%
%   OUTPUT FILENAMES: Cell array containing the filenames of all files
%                     with the given extension in the directory tree.
%
% Copyright (c) 2012 Elina Vladimirou
% Copyright (c) 2013 Jonathan W. Armond

if nargin<3
  incHidden = 0;
end
if nargin<4
  showprog = 0;
end

if showprog
  h = waitbar(0,'Finding files...');
end
filenames = {};
files = dir(directory);
for i=1:length(files)
  currentName = files(i).name;
  if ~files(i).isdir && ~isempty(regexp(currentName,regexpr, 'once')) ...
      && (incHidden || currentName(1) ~= '.')
    % A file with the prefix.
    filenames = [filenames; fullfile(directory,currentName)];
  end

  % Search subdirectories.
  if files(i).isdir && (incHidden || currentName(1) ~= '.') ...
        && ~strncmp(currentName, '..', 2)
    filenames = [filenames; kitFindFiles(fullfile(directory,currentName),regexpr)];
  end

  if showprog
    waitbar(i/length(files),h);
  end
end

if showprog
  close(h)
end

