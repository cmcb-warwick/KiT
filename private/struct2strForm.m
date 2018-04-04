function data = struct2strForm(structure)
% GETSTRUCTUREDATA Gives a cellular structure of all field names and the
% data within for use when compiling multiple datasets.
%
%    DATA = GETSTRUCTUREDATA(STRUCTURE) Produces nx2 cell structure, DATA,
%    containing the names of all sub-fields within STRUCTURE in the
%    first column, and the data found within this sub-field in the second.
%
%
% Copyright (c) 2016 C. A. Smith


% if provided data is not a structure, then output variable name and data
% therein alone
if ~isstruct(structure)
    data = who;
    data{1,2} = structure;
    return
end

% get full list of fields and subfields, then get full length
data = getFieldNames(structure);
numFields = length(data);
% get data from each of subfield
for iField = 1:numFields
    data{iField,2} = eval(['structure.' data{iField,1}]);
end



function fields = getFieldNames(structure)

% get field names
fields = fieldnames(structure);
iField = 1;
while iField < length(fields)+1
    
    if eval(['isstruct(structure.' fields{iField,1} ')'])
        
        % find number of subfields and their names
        subfields = eval(['getFieldNames(structure.' fields{iField,1} ')']);
        numSubFields = length(subfields);
        
        for iName = length(fields):-1:iField+1
            % include new names in field name list
            fields{iName-1+numSubFields,1} = fields{iName,1};
        end
        % add the field name to the start of each subfield
        homeFieldName = fields{iField,1};
        for iSub = 1:numSubFields
            subfields{iSub,1} = [homeFieldName '.' subfields{iSub,1}];
            fields{iField-1+iSub,1} = subfields{iSub,1};
        end
        
    else
        
        % if variable is not a structure, move to next field
        iField = iField+1;
    end
end
        
        


