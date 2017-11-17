function analysis=cuplMsd(analysis)
% CUPLMSD Calculate mean squared displacements
%
%   ANALYSIS = CUPLMSD(ANALYSIS) Calculates mean squared displacements for the
%   files specified in ANALYSIS. Returns same structure with results appended.
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Allocate storage.
emsd = zeros(an.maxTrackLength,3);
tmsd = zeros(an.maxTrackLength,3);

% Start point.
start = reshape(an.sisterCentreCoords(1,:), 3, an.nSisters)';

% First is always zero by definition.
for i=2:an.maxTrackLength
  % Ensemble-average MSD.
  pos = reshape(an.sisterCentreCoords(i,:), 3, an.nSisters)';
  displ = sum((pos-start).^2,2);
  emsd(i,:) = [nanmean(displ) nanstd(displ) nanserr(displ)];

  % Time-average MSD.
  displ = an.sisterCentreCoords(i:end,:) - ...
          an.sisterCentreCoords(1:end-i+1,:);
  % Break sisters into separate dimension.
  displ = reshape(displ, size(displ,1), 3, an.nSisters);
  % Compute displacements.
  displ = squeeze(sum(displ.^2,2));
  % Average across timepoint pairs.
  displ = nanmean(displ,1);
  tmsd(i,:) = [nanmean(displ) nanstd(displ) nanserr(displ)];
end

% Store result.
an.msd.ensemble.m_d = emsd(:,1);
an.msd.ensemble.s_d = emsd(:,2);
an.msd.ensemble.e_d = emsd(:,3);

an.msd.time.m_d = tmsd(:,1);
an.msd.time.s_d = tmsd(:,2);
an.msd.time.e_d = tmsd(:,3);

% Record stage.
an.stages = union(an.stages,'msd');

% Unalias analysis.
analysis = an;
