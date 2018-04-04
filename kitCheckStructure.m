function type = kitCheckStructure(struct)
% Checks the type of data structure for a given STRUCT.
%
% Copyright (c) 2017 C. A. Smith

% Initiate type.
type = 'unknown';

% Empty data structure.
if nargin<1 || isempty(struct)
  type = 'empty';
end

% Reduce any cell structure to bare bones.
while iscell(struct)
  struct = struct{1};
end

% Not a data structure.
if ~isstruct(struct)
  type = 'notstruct';
end

% Jobset or kittracking file.
if all(isfield(struct,{'kit','ROI','options'}))
  if any(isfield(struct,{'dataStruct','output'}))
    type = 'kittracking';
  else
    type = 'jobset';
  end
end

% Intra-measurements data structure.
if all(isfield(struct,{'dublVersion','microscope','plate'}))
  type = 'intrameas';
end

% Intensity data structure.
if all(isfield(struct,{'raw','norm','bg'}))
  type = 'intensity';
end

% Spot selection data structure.
if all(isfield(struct,{'selection','dataType'}))
  type = 'spotsel';
end

end