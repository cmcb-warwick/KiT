function analysis=cuplAutocorrelation(analysis)
%CUPLAUTOCORRELATION  Calculate autocorrelations
%
%   ANALYSIS = CUPLAUTOCORRELATION(ANALYSIS) Calculates autocorrelations for the
%   files specifide in ANALYSIS. Returns same structure with results appended.
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Lag times.
numLags = floor((an.maxTrackLength-1)/2);
autocorrs.t = an.time(1:numLags+1);
autocorrs.numLags = numLags;

% Compute change in sister centre position.
delX = diff(an.sisterCentreCoords(:,1:3:end));

% Compute sister breathing.
sisterSqVecs = (an.sisterCoords1 - an.sisterCoords2).^2;
sisterSep = sqrt(sisterSqVecs(:,1:3:end) + sisterSqVecs(:,2:3:end) + ...
                 sisterSqVecs(:,3:3:end));
delD = diff(sisterSep);

% Compute autocorrelation of sister centres.
autocorrs.sisters.dx = autocorrelation(delX,numLags);
autocorrs.sisters.dd = autocorrelation(delD,numLags);

% Compute autocorrelation of projected centre displacements.
autocorrs.proj.centre.dr = autocorrelation(an.proj.centre,numLags);

% Aggregate autocorrelation.
autocorrs.sisters = cuplAggregate(an,autocorrs.sisters,'dx');
autocorrs.sisters = cuplAggregate(an,autocorrs.sisters,'dd');
autocorrs.proj.centre = cuplAggregate(an,autocorrs.proj.centre,'dr');

if an.options.doTracks > 0
  % Compute change in track position.
  delX = diff(an.trackCoords(:,1:3:end));

  % Compute autocorrelation of tracks.
  autocorrs.tracks.dx = autocorrelation(delX,numLags);

  % Aggregate autocorrelation.
  autocorrs.tracks = cuplAggregate(an,autocorrs.tracks,'dx');
end

% Store result.
an.autocorrs = autocorrs;

% Record stage.
an.stages = union(an.stages,'autocorrs');

% Unalias analysis.
analysis = an;
