function pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree)
% MMFDISTPV3D calculates the p-values of inter-particle distances using a t-test
%
% SYNOPSIS pValue = mmfDistPV(maximaPos,psfSigma,numMaxima,numDegFree)
%
% INPUT  maximaPos  : Particle positions.
%        psfSigma   : Standard deviation of PSF estimate
%        numMaxima  : Number of particles.
%        numDegFree : Number of degrees of freedom for t-test.
%
% OUTPUT p-value    : Matrix of inter-particle distance p-values.
%
% REMARKS This function in its current format is not really written for
% general use but as a function to be called by mixtureModelFit.m.
%
% Copyright (c) 2011 Khuloud Jaqaman
% Copyright (c) 2012 Ed Harry
% Copyright (c) 2014 Jonathan Armond

%reserve memory for output variable
pValue = zeros(numMaxima);

%go over all pairs of particles
if size(maximaPos,2)==3
  for k = 1 : numMaxima-1
    for j = k+1 : numMaxima
      %calculate distance between the 2 particles
      x1_x2 = (maximaPos(j,1) - maximaPos(k,1));
      y1_y2 = (maximaPos(j,2) - maximaPos(k,2));
      z1_z2 = (maximaPos(j,3) - maximaPos(k,3));
      distance = sqrt(x1_x2^2+y1_y2^2+z1_z2^2);

      %get the standard deviation in the distance
      j1 = 4*(j-1)+1;
      k1 = 4*(k-1)+1;

      stdDist = x1_x2^2*(varCovMat(j1,j1) + varCovMat(k1,k1) - 2*varCovMat(j1,k1)) + ....
                y1_y2^2*(varCovMat(j1+1,j1+1) + varCovMat(k1+1,k1+1) - 2*varCovMat(j1+1,k1+1)) + ...
                z1_z2^2*(varCovMat(j1+2,j1+2) + varCovMat(k1+2,k1+2) - 2*varCovMat(j1+2,k1+2)) + ...
                2*x1_x2*y1_y2*(varCovMat(j1,j1+1) - varCovMat(j1,k1+1) - varCovMat(j1+1,k1) + varCovMat(k1,k1+1)) + ...
                2*x1_x2*z1_z2*(varCovMat(j1,j1+2) - varCovMat(j1,k1+2) - varCovMat(j1+2,k1) + varCovMat(k1,k1+2)) + ...
                2*y1_y2*z1_z2*(varCovMat(j1+1,j1+2) - varCovMat(j1+1,k1+2) - varCovMat(j1+2,k1+1) + varCovMat(k1+1,k1+2));
      stdDist = sqrt(stdDist)/distance;

      %1-sided t-test: H0: T=0, H1: T>0
      %calculate test statistic (t-distributed)
      testStat = distance/stdDist;

      %get p-value. standard deviation is 1.
      if isreal(testStat)
        pValue(j,k) = 1-tcdf(testStat,numDegFree);
      else
        pValue(j,k) = 1; % if error in varCovMat discard.
      end
    end
  end
else
  for k = 1 : numMaxima-1
    for j = k+1 : numMaxima
      %calculate distance between the 2 particles
      x1_x2 = (maximaPos(j,1) - maximaPos(k,1));
      y1_y2 = (maximaPos(j,2) - maximaPos(k,2));
      distance = sqrt(x1_x2^2+y1_y2^2);

      %get the standard deviation in the distance
      j1 = 3*(j-1)+1;
      k1 = 3*(k-1)+1;

      stdDist = x1_x2^2*(varCovMat(j1,j1) + ...
                         varCovMat(k1,k1) - 2*varCovMat(j1,k1)) ...
                + y1_y2^2*(varCovMat(j1+1,j1+1) + ...
                           varCovMat(k1+1,k1+1) - 2*varCovMat(j1+1,k1+1)) ...
                + 2*x1_x2*y1_y2*(varCovMat(j1,j1+1) - ...
                                 varCovMat(j1,k1+1) - varCovMat(j1+1,k1) + ...
                                 varCovMat(k1,k1+1));
      stdDist = sqrt(stdDist)/distance;

      %1-sided t-test: H0: T=0, H1: T>0
      %calculate test statistic (t-distributed)
      testStat = distance/stdDist;

      %get p-value. standard deviation is 1.
      if isreal(testStat)
        pValue(j,k) = 1-tcdf(testStat,numDegFree);
      else
        pValue(j,k) = 1; % if error in varCovMat discard.
      end
    end
  end
end
