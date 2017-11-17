function [outlierIdx,inlierIdx] = detectOutliers(observations,k)
%DETECTOUTLIERS detects outliers in a dataset
%
%SYNPOSIS [outlierIdx,inlierIdx] = detectOutliers(observations,k)
%
%INPUT  observations: Vector of observations (dataset).
%       k           : Roughly, for a certain value of k, observations
%                     that are k*sigma away from the mean will be
%                     considered outliers. 
%                     Optional. Default: 3.
%OUTPUT outlierIdx  : Index of observations that are considered outliers.
%       inlierIdx   : Index of observations that are considered inliers.
%
%REMARKS See Danuser 1992 or Rousseeuw & Leroy 1987 for details of
%        algorithm.
%
%Khuloud Jaqaman, October 2007

if nargin < 2 || isempty(k)
    k = 3;
end

%% outlier detection

%calculate median of observations
medObs = median(observations(:));

%get residuals, i.e. distance of observations from median
residuals = observations(:) - medObs;

%square the residuals
res2 = residuals .^ 2;

%calculate the median of the squared residuals
medRes2 = median(res2);

%define parameter to remove outliers
magicNumber2 = 1.4826^2;

%calculate test-statistic values
testValue = res2 / (magicNumber2 * medRes2);

%determine which observations are inliers and which are outliers
inlierIdx = find(testValue <= k^2);
outlierIdx = find(testValue > k^2);

%% ~~~ the end ~~~