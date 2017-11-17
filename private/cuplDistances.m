function analysis=cuplDistances(analysis)
% CUPLDISTANCES Calculate distances
%
%   ANALYSIS = CUPLDISTANCES(ANALYSIS) Calculates distances for the files
%   specified in ANALYSIS. Returns same structure with results appended.
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Position of sister tracks.
pos1 = []; pos2 = [];
for iSis = 1:an.nSisters
    pos1 = [pos1; an.sisterCoords1( : , 3*(iSis-1)+1:3*(iSis-1)+3 )];
    pos2 = [pos2; an.sisterCoords2( : , 3*(iSis-1)+1:3*(iSis-1)+3 )];
end

% Compute distance between sisters.
distances.sisters.sister.d = eudist(pos1,pos2);

% Compute distance between pair centres.
distances.sisters.pair.d = eudist(...
  an.sisterCentrePos(an.pairIdx.sisters.pair(:,1),:),...
  an.sisterCentrePos(an.pairIdx.sisters.pair(:,2),:));

% Compute distance between individual sisters.
pos = [pos1; pos2];
distances.sisters.ind.d = eudist(...
  pos(an.pairIdx.sisters.ind(:,5),:),...
  pos(an.pairIdx.sisters.ind(:,6),:));

% Aggregate distances overall.
distances.sisters.ind.m_d = nanmean(distances.sisters.ind.d,2);
distances.sisters.ind.s_d = nanstd(distances.sisters.ind.d,0,2);
distances.sisters.ind.e_d = nanserr(distances.sisters.ind.d,2);
distances.sisters.pair.m_d = nanmean(distances.sisters.pair.d,2);
distances.sisters.pair.s_d = nanstd(distances.sisters.pair.d,0,2);
distances.sisters.pair.e_d = nanserr(distances.sisters.pair.d,2);
distances.sisters.sister.m_d = nanmean(distances.sisters.sister.d,2);
distances.sisters.sister.s_d = nanstd(distances.sisters.sister.d,0,2);
distances.sisters.sister.e_d = nanserr(distances.sisters.sister.d,2);

% TODO cell means

% Store result.
an.distances = distances;

% Record stage.
an.stages = union(an.stages,'distances');

% Unalias analysis.
analysis = an;
