function [listOfNorms,normedVectors]=normList(vectors)
%calculates the norm of a list of vectors
%
%SYNOPSIS [listOfNorms,normedVectors]=normList(vectors)
%
%INPUT list of vectors (nVectorsXdimension)
%
%OUTPUT listOfNorms: list (nX1) containing the norms of the vectors
%       normedVectors: list (nXdim) containing the normed vectors
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVectors = size(vectors,1);
nDims = size(vectors,2);

%listOfNorms=zeros(nVectors,1);
listOfNorms=sqrt(sum(vectors.^2,2));

% the unit vector of length 0 is [0 0 0]
if nargout > 1
    normedVectors=zeros(size(vectors));
    normedVectors = bsxfun(@rdivide,vectors,listOfNorms);
    normedVectors(isnan(normedVectors)&~isnan(vectors)) = 0;
end

