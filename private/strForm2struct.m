function structure = strForm2struct(strForm)
% GETSTRUCTUREDATA Gives a cellular structure of all field names and the
% data within for use when compiling multiple datasets.
%
%    DATA = GETSTRUCTUREDATA(STRUCTURE) Produces nx2 cell structure, DATA,
%    containing the names of all sub-fields within STRUCTURE in the
%    first column, and the data found within this sub-field in the second.
%
%
% Copyright (c) 2016 C. A. Smith


% if provided data is not in the correct structure, give error msg
if ~iscell(strForm) || size(strForm,2) ~=2
    error('Data information provided not in the correct format.')
end

% get number of fields
numFields = size(strForm,1);
% give the data to each field
for iField = 1:numFields
    eval(['structure.' strForm{iField,1} '= strForm{iField,2};']);
end


