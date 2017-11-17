function d = meanMinDiff(A,B)
% Compute mean minimum difference between point clouds
%
% Copyright 2015 J. W. Armond

% Compute minimum difference between each point in A and set B.
sumMinDiffA = 0;
for j=1:size(A,1)
  sumMinDiffA = sumMinDiffA + min(createDistanceMatrix(A(j,:),B));
end

% Compute minimum difference between each point in B and set A.
sumMinDiffB = 0;
for j=1:size(B,1)
  sumMinDiffB = sumMinDiffB + min(createDistanceMatrix(B(j,:),A));
end

% Combine metrics to estimate difference in point clouds.
d = (sumMinDiffA + sumMinDiffB)/(size(A,1)+size(B,1));
