function [solution,resnorm,residuals,jac]=fitNGaussiansFitFun(options,x0,lb,ub,img,pixels,psfSigma)

if length(psfSigma)==3
  f = @fitNGaussians3D;
else
  f = @fitNGaussians2D;
  psfSigma = psfSigma(1);
end

[solution,resnorm,residuals,exitFlag,~,~,jac] = lsqnonlin(...
  f,x0,lb,ub,options,img,pixels,psfSigma);
switch exitFlag
    case 0
        warning('fitNGaussiansFitFun: Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals.');
    case -2
        warning('fitNGaussiansFitFun: Problem is infeasible: the bounds lb and ub are inconsistent.');
    case -4
        warning('fitNGaussiansFitFun: Line search could not sufficiently decrease the residual along the current search direction.');
end
