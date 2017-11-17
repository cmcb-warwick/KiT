function analysis=cuplSpindleLength(analysis)
% CUPLSPINDLELENGTH  Calculate spindle lengths
%
%   ANALYSIS = CUPLSPINDLELENGTH(ANALYSIS) Calculates spindle lengths for the files specified in
%   ANALYSIS. Returns same structure with results appended.
%
% Copyright (c) 2015 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Loop over cells.
cellIdx = unique(an.trackCellIdx);
nVars = 3;
d = zeros(an.nCells,nVars); % pole-pole dist, mean left pole-kt dist, mean right pole-kt dist,
for i=1:an.nCells
  % Identify poles by mean position distance from plate.
  m = nanmean(an.trackCoords(:,indexToCols(find(an.trackCellIdx==cellIdx(i)))));
  mx = m(:,1:3:end);
  poleCands = find(abs(mx)>an.options.poleCutoff);

  d(i,:) = nan(1,nVars);
  if nnz(poleCands) > 1
    [~,neg] = min(mx(poleCands));
    [~,pos] = max(mx(poleCands));
    neg = poleCands(neg);
    pos = poleCands(pos);
    negPole = m(indexToCols(neg));
    posPole = m(indexToCols(pos));
    if sign(negPole(1))~=sign(posPole(1))
      d(i,1) = eudist(negPole,posPole);

      % Extract this cells sisters.
      idx = indexToCols(find(an.sisterCellIdx==cellIdx(i)));
      c1 = nanmean(an.sisterCoords1(:,idx));
      c2 = nanmean(an.sisterCoords2(:,idx));

      % Swap sisters so x1 nearest to neg pole, x2 nearest to pos pole.
      x1 = c1(:,1:3:end);
      x2 = c2(:,1:3:end);
      idx = indexToCols(find(x1>x2));
      [c2(:,idx), c1(:,idx)] = deal(c1(:,idx), c2(:,idx));

      % Compute mean distances.
      d(i,2) = nanmean(eudist(negPole,reshape(c1,[],3)));
      d(i,3) = nanmean(eudist(posPole,reshape(c2,[],3)));
    end
  end
end

% Store result.
an.spindle.length = d(:,1);
an.spindle.leftKfibre = d(:,2);
an.spindle.rightKfibre = d(:,3);

% Record stage.
an.stages = union(an.stages,'spindle');

% Unalias analysis.
analysis = an;
