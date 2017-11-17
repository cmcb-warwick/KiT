function [c,distance] = clusterSpot(coordList,cutoff,correction)
% CLUSTERSPOT cluster spots based on the theoretical or corrected psf size
%
% SYNOPSIS: cluster = clusterSpot(coordList,cutoff,correction)
%
% INPUT coordList : [nSpot,nDim] coordinate list of the spots to be cluster.
%		cutoff : Distance separating 2 clusters.
%		correction(opt) : 1 - Calculate the correction of the psf
%                                 0 - Use directly the psf size as pass in argument.
%                         def : 0
%
% OUTPUT c : 1-nSpot Vector with the cluster index for each spot coordinate.
%        distance : 1-by-nSpot-1 Vector with the distance of the clusters
% REMARKS
%
% SEE ALSO linkage,cluster
%
% EXAMPLE
%         cl = clusterSpot(coordList,sigma);
%
% created by: Jacques Boisvert DATE: 24-May-2012
% Revision by EHarry 19/11/12, if working in 3D will warp z coordinates so can
% use a single distance cutoff (assuming psf in x and y are the same)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ARGUMENTS VALIDATION
if nargin < 1 || isempty(coordList)
    error('coordList is missing');
end

if nargin < 2
    error('psf size is missing');
end

if nargin < 3
    %Default
    correction = 0;
end
nSpot = size(coordList,1);
%------------------------------------------------------------------------

%% PSF CORRECTION
if correction
end

%% COORDINATE WARPING IN Z
if length(cutoff) == 3
    coordList(:,3) = coordList(:,3) .* cutoff(1) ./ cutoff(3);
    cutoff = cutoff(1);
end

%% HIERERCHICAL CLUSTERING

if nSpot >1
    % Regroup cluster based on cutoff.
    z = linkage(coordList,'single');
    c = cluster(z,'cutoff',cutoff,'criterion','distance');
    distance = z(:,3);
else
    c = 1;
    distance = NaN;
end
%------------------------------------------------------------------------
