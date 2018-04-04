function data = combineStrForms(strForm1,strForm2)
% GETSTRUCTUREDATA Gives a cellular structure of all field names and the
% data within for use when compiling multiple datasets.
%
%    DATA = GETSTRUCTUREDATA(DATA1,DATA2) Produces nx2 cell structure, DATA,
%    containing the names of all sub-fields within STRUCTURE in the
%    first column, and the data found within this sub-field in the second.
%
%
% Copyright (c) 2016 C. A. Smith


% if provided data is not in the correct structure, give error msg
if ~iscell(strForm1) || size(strForm1,2) ~=2
    error('Data information provided in first entry not in the correct format.')
elseif ~iscell(strForm2) || size(strForm2,2) ~=2
    error('Data information provided in second entry not in the correct format.')
end

% get list of structure components
components = strForm1(:,1);
numComps = size(strForm1,1);
% push data1 to output data variable
data = strForm1;

% loop through original components and concatenate
for iComp = 1:numComps
    jComp = find(strcmp(strForm2(:,1),components(iComp)),1);
    if ~isempty(jComp)
        data{iComp,2} = [data{iComp,2}; strForm2{jComp,2}];
    end
end

% loop through new components and alert user of any not present in original
components = strForm2(:,1);
numComps = size(strForm2,1);
for jComp = 1:numComps
    iComp = find(strcmp(strForm1(:,1),components(jComp)),1);
    if isempty(iComp)
        warning('Structural component, %s, is present only in new data.',components{jComp});
    end
end
