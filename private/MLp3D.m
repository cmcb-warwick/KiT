function [paramsf, conflim] = MLp3D(data, params_start,verbose,numbins)
% MLp3D maximum likelihood fit to p3D(r)
%   [PARAMSF, CONFLIM] = MLp3D(DATA, PARAMS_START, NUMBINS) takes the vector, 
%   DATA, of Euclidean distances calculated in a 3D geometry and the  
%   best-guess set of parameters, PARAMS_START, and performs a maximum 
%   likelihood fit to DATA with the theoretical probability distribution, p3D(r):
%        p3D(r) = ((2/pi))^.5*exp(-(x.^2+mu^2)./(2*sigma^2)).*x.*sinh(mu*x./sigma^2)./(sigma*mu). 
%   A histogram of the binned data set
%   into NUMBINS bins along with the fit curve is shown. PARAMSF gives the 
%   fit distance as the first element and the fit RMSD as the second element.
%   The 68% confidence limits of the fit distance estimator is given in CONFLIM. 
%   The first element is the error to the left of the found distance, and 
%   the second element is the error to the right of the found distance. A 
%   plot of the log likelihood for the distance estimator is plotted
%   and should appear approximately parabolic around the peak for the 
%   calculated errors to have merit. 
%
%   Written by L. Stirling Churchman, Stanford University, 2005.
%   Version 1.0
%    Modified by E.D. Salmon, Oct2016 as described above.
%   Reference:
%     Barlow RJ. 1989. Statistics: a guide to the use of statistical
%     methods in the physical sciences. New York, NY: John Wiley and Sons.
%     204 p.
%   Edits:
%     Here Eqn (6) is substituted for Eqn 4 in Churchman et al Biophysics J
%     2006 for 3D fit. I did not change names of functions; just equation
%     in last 2 functions. Ted Salmon. October, 2016
%     All NaNs now removed from data prior to processing. Option to
%     suppress plotting and display of results. Chris Smith. January, 2018

if nargin<3 || isempty(verbose)
  verbose = 0;
end
if nargin<4
  numbins = 25;
end

% can't handle nans or zeros
data( isnan(data) | data==0 ) = [];
num_pts = length(data);
 
% Define the likelihood function and maximize it
fun = @(params)(-prod(100*P3D(params, data)));
options = optimset('TolFun', 1e-5);
%options = optimset(options,'Display','iter')
paramsf = fminsearch(fun, params_start, options);
 
% View the results of the maximum likelihood fit
if verbose
    figure();
    subplot(2,1,1);
    [N, xbin] = hist(data, numbins);
    b = bar(xbin, N/trapz(N),'FaceColor',[240 248 255]/255);
    hold on
    [~, xbin] = hist(data, 50);
    deltx = xbin(6)-xbin(5);
    ex_pts = floor((xbin(1)-deltx/2)/deltx);
    pts = (xbin(1)-ex_pts*deltx):deltx:(xbin(1)-deltx);
    Q=intP3D(paramsf, ((xbin(1)-deltx/2)-ex_pts*deltx), xbin(end)+deltx/2, deltx);
    plot([pts xbin], 50/numbins*Q, 'r-','LineWidth',2);
    line([paramsf(1) paramsf(1)],ylim,'Color',[0 0 1],'LineWidth',2);
    xlabel('Distance measurements')
    ylabel('Frequency')
    title('Histogram of data with the fit to \itp_3_D(r)')
    xlim([0 200]);
    set(gca,'FontSize',16);
    muleg = sprintf('mu = %.2f nm',paramsf(1));
    legend('Data', 'ML fit to \itp_3_D(r)',muleg);
    legend boxoff
    hold off
end
 
% Determine the confidence limits of the fit distance parameter
 sigma_range = 4*paramsf(2)/sqrt(num_pts);
resolution = 0.001*paramsf(2);
if (paramsf(2)<paramsf(1))
    x= ((paramsf(1)-sigma_range):resolution:(paramsf(1)+sigma_range));
else
    x= (resolution:resolution:(paramsf(1)+sigma_range));
end
% Calculate the log likelihood for a range of distances    
x_pts = length(x);
L = zeros(x_pts, 1);
for i = 1:x_pts
    L(i)=sum(log(P3D([x(i) paramsf(2)],data)));
end
 
% Plot the log likelihood
if verbose
    subplot(2,1,2);
    plot(x, L,'b','LineWidth',2);
    grid on
    xlabel('Distance')
    ylabel('Log Likelihood')
    title('Log Likelihood versus the distance estimator')
    set(gca,'FontSize',16);
end
% Determine the distances one sigma from the center distance as the 68%
% confidence limits
Ls = sort(abs(L-(max(L)-0.5)));
f1 = find(abs(L-(max(L)-0.5))==Ls(1));
f2 = find(abs(L-(max(L)-0.5))==Ls(2));
if (f1 < f2)
    conflim(1) = paramsf(1)-x(f1);
    conflim(2) = x(f2) - paramsf(1);
else
    conflim(1) = paramsf(1)-x(f2);
    conflim(2) = x(f1) - paramsf(1);
end
 
if verbose
  % Display the results
  disp('Results of the maximum likelihood fit to p3D(r):');
  disp(['   n = ' num2str(num_pts)]);
  disp(['   Distance ' num2str(paramsf(1)) ' + ' num2str(conflim(2)) ' - ' num2str(conflim(1))]);
  disp(['   Sigma ' num2str(paramsf(2))]);
end
 
end

%% Sub-functions

function output = P3D(params, x)
% P3D Codes for p3D(r) from Eq. (6)
%   OUTPUT = P3D(PARAMS, X)For the given set of PARAMS and a list of x
%   values, the p3D(x) is found.
mu = params(1);
sigma = params(2);
output = ((2/pi))^.5*exp(-(x.^2+mu^2)./(2*sigma^2)).*x.*sinh(mu*x./sigma^2)./(sigma*mu);

end
 
 
function Q = intP3D(params, x1, x2, deltx)
% INTP3D Integrated p3D(r)
%   Q = INTP3D(PARAMS, X1, X2, DELTX) integrates the probability
%   described by p3D(r) through steps of DELTAX from X1 to X2.
 
num = round((x2-x1)/deltx);
mu = params(1);
sigma = params(2);
F = @(x) ((2/pi)^.5)*exp(-(x.^2+mu^2)./(2*sigma^2)).*x.*sinh(mu*x./sigma^2)./(sigma*mu);
for i = 1:num
    Q(i) = integral(F, x1+(i-1)*deltx,x1+i*deltx);
end

end
