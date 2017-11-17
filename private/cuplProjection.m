function analysis=cuplProjection(analysis)
% CUPLPROJECTION Run analyses on selected files
%
%   ANALYSIS = CUPLPROJECTION(ANALYSIS) Compute projection of displacement
%   vector onto sister-sister vector at previous timepoint. Potentially
%   useful if no coordinate-system available.
%
% Copyright (c) 2013 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Sister-sister vectors.
sisVecs = an.sisterCoords2(1:end-1,:) - an.sisterCoords1(1:end-1,:);
% Normalize.
for i=1:an.nSisters
  cols = indexToCols(i);
  sisVecs(:,cols) = bsxfun(@rdivide,sisVecs(:,cols),...
                           sqrt(sum(sisVecs(:,cols).^2,2)));
end

% Sister centre displacement vectors.
centreDisp = diff(an.sisterCentreCoords);

% Projected displaced is dot product of centre disp vector and normalized
% sister-sister vector.
projCentreDisp = zeros(an.maxTrackLength-1,an.nSisters);
for i=1:an.nSisters
  cols = indexToCols(i);
  projCentreDisp(:,i) = dot(centreDisp(:,cols),sisVecs(:,cols),2);
end

% Store results.
an.proj.centre = projCentreDisp;

% Record stage.
an.stages = union(an.stages,'projection');

% Unalias analysis.
analysis = an;
