function a=structCopyMissingFields(a,b)
% STRUCTCOPYMISSINGFIELDS Copy fields that B has but A does not to A

if ~(isstruct(a) && isstruct(b))
  error('A and B must be structs');
end

aFields = fieldnames(a);
bFields = fieldnames(b);

% Recursive if any field of A is a struct.
for i=1:length(aFields)
  f = aFields{i};
  if isstruct(a.(f)) && ismember(f,bFields);
    a.(f) = structCopyMissingFields(a.(f), b.(f));
  end
end

% Missing fields.
missing = setdiff(bFields,aFields);

% Copy over.
for i=1:length(missing)
  f = missing{i};
  a.(f) = b.(f);
end
