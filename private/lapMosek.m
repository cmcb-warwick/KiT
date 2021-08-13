function [x, y] = lapMosek(cc, NONLINK_MARKER, augmentCC, noLinkCost)
%
%LAPMOSEK solves the linear assignment problem for a given cost matrix
%using the mosek optimization library, which requires separate installation
%See: https://www.mosek.com/ Free academic licenses available 
%
% A linear assignment tries to establish links between points in two sets.
% One point in set A can only link to one point in set B and vice versa.
% The cost associated with the link from element i of A to element j of B
% is given by cc(i,j).
%
% SYNOPSIS [x, y] = lapMosek(cc, NONLINK_MARKER, augmentCC, noLinkCost)
%
% INPUT:  cc: cost matrix, which has to be square. Set cost(i,j) to the
%             value of NONLINK_MARKER, if the link is not allowed. cc can
%             be input as a sparse matrix (in which case there won't be any
%             NONLINK_MARKERs).
%             Elements of the cost matrix cannot be inf or nan!
%
%             For an unequal number of points (or, generally, if birth and
%             death is to be allowed) the cost matrix is formed as
%             a 2-by-2 catenation of four sub-matrices. If there are n
%             elements in set A and m elements in set B, sub-matrix (1,1)
%             is an n-by-m matrix with the cost for each link. The two
%             off-diagonal sub-matrices make non-links possible. They are
%             both diagonal square matrices (2,1: m-by-m; 1,2: n-by-n) with
%             the cost for not linking a point (e.g. determined by a
%             maximum search radius) in the diagonal and NONLINK_MARKER off
%             the diagonal. Sub-matrix (2,2) is a m-b-n matrix of any
%             not-NONLINK_MARKER value (suggested: one, except for sparse
%             matrix input)
%
% NONLINK_MARKER : value to indicate that two points cannot be linked.
%             Default: -1. NaN is not allowed here.
%
% extendedTesting (used to be optional, now input has no effect and
%                  extendedTesting will always be performed)
%                LAP will always make sure that:
%                  - There cannot be NaNs in cc
%                  - In every row and every column of cc there must be
%                    at least 1 element that is not a NONLINK_MARKER.
%
% augmentCC  (optional, [{0}/1]) If 1, the cost matrix will be augmented to
%            allow births and deaths.
% noLinkCost (optional, [{maximum cost + 1}]) Cost for linking a feature
%            to nothing.
%
%
% OUTPUT: x: The point A(i) links to B(x(i))
%         y: The point B(j) links to A(y(j))
%
%            Any x > m or y > n indicates that this point is not linked.
% Jonathan U Harrison 2019-18-01 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=====================
% TEST INPUT
%=====================

if (nargin == 1) || isempty(NONLINK_MARKER)
    NONLINK_MARKER = -1;
elseif isnan(NONLINK_MARKER)
    error('NONLINK_MARKER cannot be NaN!')
end

if nargin < 4 || isempty(augmentCC)
    augmentCC = 0;
end

% test size
scc = size(cc);
% check size only if no augmentation
if ~augmentCC

    if scc(1) ~= scc(2) || length(scc) > 2
        error('cost must be a 2D square matrix!')
    end

elseif length(scc) > 2

    error('cost must be a 2D matrix!')

else

    % if we're augmenting, sparse matrices are produced. This will be
    % problematic with cost 0
    if any(cc(:)==0) && ~issparse(cc)
        validCC = cc ~= NONLINK_MARKER;
        if any(any(cc(validCC))) < 0
            warning('LAP:negativeCosts',...
                'there are negative and zero costs. This could lead to errors')
        end
        cc(validCC) = cc(validCC) + 10;
    end

end

% extended testing has been made mandatory

if nargin < 4 || isempty(noLinkCost)
    noLinkCost = max(max(cc)) + 1;
end

% if we have -1 as non-link-marker, and we get it everywhere, we get
% a segmentation fault.
if noLinkCost == NONLINK_MARKER
    noLinkCost = noLinkCost + 1;
    if noLinkCost == 0
        noLinkCost = 0.01;
    end
end

%%%%%%%%%%%%%%%%%%

%=================================
% AUGMENT COST MATRIX IF SELECTED
%=================================
if augmentCC

    % expand the m-by-n cost matrix to an (m+n)-by-(n+m) matrix, adding
    % diagonals with noLinkCost

    % in the lower right corner, we want the lowest cost of all - take
    % minimum cost, subtract 5 and make sure it's not 0.
    % NONLINK_MARKER is not a problem because we make the matrix sparse
    % before augmenting.
    %
    % not needed anymore since we use cc' as costLR
%     minCost = min(min(cc(cc~=NONLINK_MARKER)))-5;
%     if minCost == 0
%         minCost = -2;
%     end


    % check if sparse
    if issparse(cc)

        % mmDiag = spdiags(noLinkCost * ones(scc(1),1), 0, scc(1), scc(1));
        % nnDiag = spdiags(noLinkCost * ones(scc(2),1), 0, scc(2), scc(2));
        % nmMat  = sparse(ones(scc(2), scc(1)));
        % cc = [cc, mmDiag; nnDiag, nmMat];

        % use transposed costmatrix for LR - otherwise, nonLinkCost has no
        % effect
        %[rowIdx, colIdx] = find(cc);
        %costLR = sparse(colIdx,rowIdx,minCost*ones(length(colIdx),1),scc(2),scc(1));
        costLR = cc';

        cc = [cc, spdiags(noLinkCost * ones(scc(1),1), 0, scc(1), scc(1));...
            spdiags(noLinkCost * ones(scc(2),1), 0, scc(2), scc(2)),...
            costLR];

    else

        % mmDiag = diag(noLinkCost * ones(scc(1),1));
        % nnDiag = diag(noLinkCost * ones(scc(2),1));
        % nmMat  = sparse(ones(scc(2), scc(1)));
        % cc = [cc, mmDiag; nnDiag, nmMat];

        % use transposed costmatrix for LR - otherwise, nonLinkCost has no
        % effect
        %[rowIdx, colIdx] = find(cc ~= NONLINK_MARKER);
        %costLR = sparse(colIdx,rowIdx,minCost*ones(length(colIdx),1),scc(2),scc(1));

        % make cc sparse. Take NLM in cc into account!
        cc(cc==NONLINK_MARKER) = 0;
        cc=sparse(cc);
        costLR = cc';

        cc = [cc, diag(noLinkCost * ones(scc(1),1)); ...
            diag(noLinkCost * ones(scc(2),1)), ...
            costLR];

    end

    % remember that the size of the matrix has increased!
    scc = [sum(scc), sum(scc)];

    clear costLR colIdx rowIdx

end
%%%%%%%%%%%
cc = cc';
% find the significant elements. If sparse input, find nonzero elements
if issparse(cc)
    [rowIdx, colIdx, val] = find(cc);
else
    [rowIdx, colIdx] = find(cc ~= NONLINK_MARKER);
    val = cc(cc ~= NONLINK_MARKER);
end

nonLinkIdx = setdiff(1:numel(cc),sub2ind(scc,rowIdx,colIdx));
[nlmRow, nlmCol] = ind2sub(scc, nonLinkIdx);

% test that all cols and all rows are filled, and that there are no nans
    if issparse(cc)
        if (~all(sum(cc,1)) || ~all(sum(cc,2)))
            error('LAP:BadCostMatrix',...
                'Rows and columns of the cost matrix must allow at least one possible link!')
        end
    else
        if (~all(sum(cc~=NONLINK_MARKER,1)) || ~all(sum(cc~=NONLINK_MARKER,2)))
            error('LAP:BadCostMatrix',...
                'Rows and columns of the cost matrix must allow at least one possible link!')
        end
    end

    if any(~isfinite(val))
        error('LAP:NanCostMatrix','Cost matrix cannot contain NaNs or Inf!')
    end

%%%%%%%%%%%%%%
import mosek.fusion.*;

M = Model('lapMosek');  
xx = M.variable([scc(1),scc(2)], Domain.binary());

%each value only once per dim
for d = 1:2
    M.constraint(Expr.sum(xx,d-1), Domain.lessThan(2.));
    M.constraint(Expr.sum(xx,d-1), Domain.greaterThan(1.));
end
%add constraint for non link markers
for i = 1:length(nlmRow)
        M.constraint(xx.index(nlmRow(i),nlmCol(i)),Domain.equalsTo(0.0));
end
% Set max solution time
M.setSolverParam('mioMaxTime', 60.0);
% Set max relative gap (to its default value)
M.setSolverParam('mioTolRelGap', 1e-4);
% Set max absolute gap (to its default value)
M.setSolverParam('mioTolAbsGap', 0.0);

% Set the objective function to (c^T * x)
M.objective('obj', ObjectiveSense.Minimize, Expr.dot(cc, xx));

% Solve the problem
M.solve();
% Get the solution values
sol = xx.level();

[x,y,~] = find(reshape(round(sol),scc)); %output indices of chosen connections
x = int32(x);
y = int32(y);
% totalcost = 0;
% for i=1:scc(1)/2
%    totalcost = totalcost + cc(x(i),i);
% end
% totalcost
