function data=cuplAggregate(an,data,field)
% CUPLAGGREGATE Aggregates data and stores as new fields in data
%
% Copyright (c) 2013 Jonathan Armond

data.(['m_' field]) = nanmean(data.(field),2);
data.(['s_' field]) = nanstd(data.(field),0,2);
data.(['e_' field]) = nanserr(data.(field),2);

data.(['cm_' field]) = cuplCellMean(an,data.(field));
data.(['cs_' field]) = cuplCellStd(an,data.(field));
data.(['ce_' field]) = cuplCellSerr(an,data.(field));
