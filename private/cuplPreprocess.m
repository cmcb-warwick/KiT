function analysis=cuplPreprocess(analysis)
% CUPLPREPROCESS Run analyses on selected files
%
%   ANALYSIS = CUPLPREPROCESS(ANALYSIS) Preprocesses files to count number of
%   NaNs, number of tracks per cell, etc..
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Allocate and initialise storage.
an.nSistersPerCell = zeros(an.nLoadedCells,1);
an.nTrackPerCell = zeros(an.nLoadedCells,1);

% Preprocess data.

% Count NaNs. NaNs are always mirrored between sisters and coordinates.
an.sisterNanCount = sum(isnan(an.sisterCoords1(:,1:3:end)),1);
an.trackNanCount = sum(isnan(an.trackCoords(:,1:3:end)),1);

% Sister pair is accepted if less than or equal to percentNan NaN fraction.
% Binary vector: 1 for accepted pair. These are preliminarily selected,
% prior to cell acceptance.
maxNans = an.options.percentNan*an.maxTrackLength;
prelimAcceptedSisters = (an.sisterNanCount <= maxNans);
prelimAcceptedTracks = (an.trackNanCount <= maxNans);

for i=1:an.nLoadedCells
  an.nSistersPerCell(i) = sum(prelimAcceptedSisters(an.sisterCellIdx==i));
  an.nTracksPerCell(i) = sum(prelimAcceptedTracks(an.trackCellIdx==i));
end

% Accept cell if greater than or equal to mininum number of accepted
% pairs. Combine with inspectionAcceptedCells.
% Binary vector: 1 for accepted pair
autoAcceptedCells = an.nSistersPerCell >= an.options.minSistersPerCell;
an.acceptedCells = autoAcceptedCells; % & an.inspectionAcceptedCells;
an.nSistersPerCell = an.nSistersPerCell(an.acceptedCells);
an.nTracksPerCell = an.nTracksPerCell(an.acceptedCells);
an.nCells = sum(an.acceptedCells);
an.nSisters = sum(an.nSistersPerCell);
an.nTracks = sum(an.nTracksPerCell);

% Fixup acceptedSisters index to account for only accepted cells.
an.acceptedCellsIdx = find(an.acceptedCells);
an.acceptedSisters = ismember(an.sisterCellIdx, an.acceptedCellsIdx) & ...
    prelimAcceptedSisters;
% Fixup sisterCellIdx to point to accepted cells.
an.cellIdxMap = zeros(1,an.nLoadedCells);
an.cellIdxMap(an.acceptedCells) = 1:an.nCells;
an.sisterCellIdx = an.cellIdxMap(an.sisterCellIdx(an.acceptedSisters));

% Fixup accepted... for tracks.
if an.hasTracks
  an.acceptedTracks = ismember(an.trackCellIdx, an.acceptedCellsIdx) & ...
      prelimAcceptedTracks;
  an.nAcceptedTracks = sum(an.acceptedTracks);
  % Fixup trackCellIdx to point to accepted cells.
  an.trackCellIdx = an.cellIdxMap(an.trackCellIdx(an.acceptedTracks));
end

% Remove unaccepted sisters/tracks from data.
acceptedSistersCols = indexToCols(find(an.acceptedSisters));
an.sisterCoords1 = an.sisterCoords1(:,acceptedSistersCols);
an.sisterCoords2 = an.sisterCoords2(:,acceptedSistersCols);
if an.hasTracks
  acceptedTracksCols = indexToCols(find(an.acceptedTracks));
  an.trackCoords = an.trackCoords(:,acceptedTracksCols);
end
% if an.hasTrackInt
%   for i = 1:an.nChannels
%     % FIXME Kit doesn't seem to calculate for all tracks...
%     an.trackInt(i).mean =   an.trackInt(i).mean(:,an.acceptedTracks);
%     an.trackInt(i).max =    an.trackInt(i).max(:,an.acceptedTracks);
%     an.trackInt(i).min =    an.trackInt(i).min(:,an.acceptedTracks);
%     an.trackInt(i).median = an.trackInt(i).median(:,an.acceptedTracks);
%   end
% end

% Find sister centres. Tx3N
an.sisterCentreCoords = 0.5 * (an.sisterCoords1 + an.sisterCoords2);

% Find sister centre mean position. Nx3
an.sisterCentrePos = reshape(nanmean(an.sisterCentreCoords,1),3,an.nSisters)';

if an.hasTracks
  % Find track mean position. Nx3
  an.trackPos = reshape(nanmean(an.trackCoords,1),3,an.nTracks)';
end

% Allocate storage.
an.sisterRadii = zeros(an.nSisters,1);
an.trackRadii = zeros(an.nTracks,1);
an.cellCentres = zeros(an.maxTrackLength,3*an.nCells);
an.meanCellCentres = zeros(an.nCells,3);

% Find cell centres (using all sisters).
for i=1:an.nCells
  accIdx = an.sisterCellIdx == i;
  sisterCols = indexToCols(find(accIdx));
  cellCols = indexToCols(i);
  an.cellCentres(:,cellCols) = [...
    nanmean(an.sisterCentreCoords(:,sisterCols(1:3:end)),2),... % X
    nanmean(an.sisterCentreCoords(:,sisterCols(2:3:end)),2),... % Y
    nanmean(an.sisterCentreCoords(:,sisterCols(3:3:end)),2)];   % Z

  an.meanCellCentres(i,:) = nanmean(an.cellCentres(:,cellCols));

  switch an.options.correctDrift
    case 1
      % Subtract off average centre.
      reppedCentres = repmat(nanmean(an.cellCentres(:,cellCols)),1,sum(accIdx));
      an.sisterCoords1(:,sisterCols) = ...
          bsxfun(@minus,an.sisterCoords1(:,sisterCols),reppedCentres);
      an.sisterCoords2(:,sisterCols) = ...
          bsxfun(@minus,an.sisterCoords2(:,sisterCols),reppedCentres);
      an.sisterCentreCoords(:,sisterCols) = ...
          bsxfun(@minus,an.sisterCentreCoords(:,sisterCols),reppedCentres);
    case 2
      % Subtract off centre every frame.
      reppedCentres = repmat(an.cellCentres(:,cellCols),1,sum(accIdx));
      an.sisterCoords1(:,sisterCols) = an.sisterCoords1(:,sisterCols) - reppedCentres;
      an.sisterCoords2(:,sisterCols) = an.sisterCoords2(:,sisterCols) - reppedCentres;
      an.sisterCentreCoords(:,sisterCols) = an.sisterCentreCoords(:,sisterCols) -...
          reppedCentres;
  end

  % Find sister radii.
  sisIdx = an.sisterCellIdx == i;
  an.sisterRadii(sisIdx) = eudist(an.sisterCentrePos(sisIdx,:), ...
                               repmat(an.meanCellCentres(i,:),sum(sisIdx),1));

  if an.hasTracks
    % Find track radii.
    trIdx = an.trackCellIdx == i;
    an.trackRadii(trIdx) = eudist(an.trackPos(trIdx,:), ...
                                  repmat(an.meanCellCentres(i,:),sum(trIdx),1));
  end
end

% Precompute pairing indices (using accepted tracks/sisters).
nSisterPairings = sum(arrayfun(@(n) condnchoosek(n,2), ...
                               an.nSistersPerCell));
an.pairIdx.sisters.ind = zeros(4*nSisterPairings,6,'uint32'); % uses counter 1
an.pairIdx.sisters.pair = zeros(nSisterPairings,2,'uint32');  % uses counter 2
if an.hasTracks
  nTrackPairings = sum(arrayfun(@(n) condnchoosek(n,2), ...
                               an.nTracksPerCell));
  an.pairIdx.tracks.ind = zeros(nTrackPairings,2,'uint32');   % uses counter 3
end
c = ones(3,1); % counters
for i=1:an.nCells
  sisters = find(an.sisterCellIdx == i);
  for j=1:length(sisters)
    for k=j+1:length(sisters)
      % Individual tracks of this to non-sister track.
      an.pairIdx.sisters.ind(c(1)  ,1:4) = [sisters([j k]) 1 1]; % coords1/1 tracks
      an.pairIdx.sisters.ind(c(1)+1,1:4) = [sisters([j k]) 1 2]; % coords1/2 tracks
      an.pairIdx.sisters.ind(c(1)+2,1:4) = [sisters([j k]) 2 1]; % coords2/1 tracks
      an.pairIdx.sisters.ind(c(1)+3,1:4) = [sisters([j k]) 2 2]; % coords2/2 tracks
      c(1) = c(1)+4;

      % Sister pair to sister pair.
      an.pairIdx.sisters.pair(c(2),:) = sisters([j k]);
      c(2) = c(2)+1;
    end
  end

  if an.hasTracks
    tracks = find(an.trackCellIdx == i)';
    for j=1:length(tracks)
      for k=j+1:length(tracks)
        an.pairIdx.tracks.ind(c(3),:) = [tracks([j k])];
        c(3) = c(3)+1;
      end
    end
  end
end

% Add indexes into stacked coords1/coords2 to an.pairIdx.sisters.ind.
% This is to allow computing things pairwise on e.g. x-coordinates as
%  xcoords = [an.sisterCoords1(:,1:3:end) an.sisterCoords2(:,1:3:end)];
% and indexing via an.pairIdx.sisters.ind(:,5,6);
an.pairIdx.sisters.ind(:,5:6) = an.pairIdx.sisters.ind(:,1:2) + ...
    (an.pairIdx.sisters.ind(:,3:4)-1)*an.nSisters;

% Compute displacement projections.
an = cuplProjection(an);

% Record stage.
an.stages = union(an.stages,'preprocess');

% Unalias analysis.
analysis = an;

end


function b=condnchoosek(n,k)
% Returns nchoosek if n >= k, otherwise 0.

if n<2
  b = 0;
else
  b = nchoosek(n,2);
end

end
