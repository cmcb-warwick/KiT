function trueOrFalse = isApproxEqual(A,B,alpha,relOrAbs)
%ISAPPROXEQUAL checks whether two values (matrices) are separated by less than a certain tolerance
%
% Remarks: The algorithm tests for relative difference. For each
%           element, the difference is expressed as ratio of the smaller
%           to the larger value, i.e. the relative difference is calculated
%           relative to the larger value.
%
% SYNOPSIS trueOrFalse = isApproxEqual(A,B,alpha,relOrAbs)
%
% INPUT    A,B  :  values to be compared. If either is an array, the other has
%                   to be either scalar or an array of the same size.
%          alpha: (opt) tolerance by which the two values are allowed to
%                   differ. Default: 0.001. Has to be between 0 and 1
%          relOrAbs: (opt) 'relative', if tolerance is relative, 'absolute'
%                   if absolute tolerance. Default: 'relative'
% OUTPUT   trueOrFalse : logical array of size of A,B with 1 for every true
%                           assigment
%
%c: 05/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============
% TEST INPUT
%=============

% argin
if nargin < 2
    error('please specify at least two input arguments')
end

% alpha
if nargin < 3 || isempty(alpha)
    alpha = 0.001;
else
    alpha = abs(alpha);
    if alpha < 0 || alpha > 1
        error('alpha has to be between 0 and 1')
    end
end

if nargin < 4 || isempty(relOrAbs)
    relOrAbs = 'relative';
end

% A and B
if ~isnumeric(A) || ~isnumeric(B)
    error('please specify numeric input for A and B')
end

% return false if one of the two is empty
if isempty(A) || isempty(B)
    trueOrFalse = false;
    return
end

sizA = size(A);
sizB = size(B);

scalarA = length(sizA) == 2 & all(sizA(1:2)==[1,1]);
scalarB = length(sizB) == 2 & all(sizB(1:2)==[1,1]);

switch scalarA + 2*scalarB
    % 0: none scalar
    % 1: A scalar
    % 2: B scalar
    % 3: all scalar

    case 0
        if any(sizA ~= sizB)
            error('non-scalar A and B must have the same size')
        end

    case 1

        sizA = sizB;
        A = repmat(A,sizA);

    case 2

        sizB = sizA;
        B = repmat(B,sizB);

    case 3

        % all ok
end



%==========================


%==================
% COMPARE
%==================

% turn off warnings
oldWarnings = warning;
warning off MATLAB:divideByZero

% test for both zero - otherwise, 0==0 is false
bothZero = (A==0) & (B==0);



switch relOrAbs
    case 'relative'
        % for relative comparison: relative difference smaller alpha OR both zero. The
        % min operator takes care of different signs, too, because a large number
        % will be added to one if A,B, differ in their sign
        trueOrFalse = (1-min((B./A),(A./B))  <=alpha)  | bothZero;

    case 'absolute'
        % check whether the absolute difference between A and B is smaller
        % than alpha
        trueOrFalse = abs(A-B) <= alpha;
        
    otherwise
        error('input argument relOrAbs has to be either ''relative'' or ''absolute'', not %s',relOrAbs)
end


% turn warnings back on
warning(oldWarnings);