function analysis=cuplStripData(analysis)
% CUPLSTRIPDATA Strip excessive data from stored analysis
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

if nargin<1
    error('No analysis struct supplied.');
end

% Alias analysis.
an = analysis;

% Pair indexing.
if isfield(an,'pairIdx')
  an = rmfield(an,'pairIdx');
end

% Crosscorrelation pairs.
if isfield(an,'crosscorrs')
  if isfield(an.crosscorrs.sisters.sister,'dx')
    an.crosscorrs.sisters.sister = rmfield(an.crosscorrs.sisters.sister,{'dx','x'});
  end
  if isfield(an.crosscorrs.sisters.pair,'dx')
    an.crosscorrs.sisters.pair = rmfield(an.crosscorrs.sisters.pair,{'dx','x'});
  end
  if isfield(an.crosscorrs.sisters.ind,'dx')
    an.crosscorrs.sisters.ind = rmfield(an.crosscorrs.sisters.ind,{'dx','x'});
  end
end

% Record stage.
an.stages = union(an.stages,'stripped');

% Unalias analysis.
analysis = an;
