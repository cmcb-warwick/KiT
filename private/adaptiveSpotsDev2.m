function [spots,threshold] = adaptiveSpotsDev2(movie,lambda,...
                                           realisticNumSpots, ... 
                                           flatBackground, verbose, nPieces)
% ADAPTIVESPOTS Adaptive thresholding for spot detection.
%
% Copyright 2015 J. W. Armond
%%Adaptive spots development: theshold that varies in time
%
%
%JUH
%%%%%%%%%%


if nargin<2 || isempty(lambda)
  lambda = 0;
end

if nargin<3 || isempty(realisticNumSpots)
  realisticNumSpots = 100;
end

if nargin<4 || isempty(flatBackground)
    flatBackground = 0;
end

if nargin<5 || isempty(verbose)
  verbose = 0;
end

if nargin<6 || isempty(nPieces)
  nPieces = 1;
end

have92 = ~verLessThan('images','9.2');

% Go over all frames and find local maxima.
nFrames = size(movie,4);
locs = cell(nFrames,1);
meanInt = zeros(nFrames,1);
for i=1:nFrames
  img = movie(:,:,:,i);
  meanInt(i) = mean(img(:));
  % TODO options?
  if have92
    imgF = imgaussfilt3(img,2,'FilterSize',3);
    bkgd = imgaussfilt3(img,16,'FilterSize',63);
  else
    imgF = fastGauss3D(img,2,3);
    bkgd = fastGauss3D(img,16,63);
  end
  if flatBackground
  bkgd = ones(size(imgF))*mean(bkgd(:));
  end
  amp = imgF-bkgd;

  bw = imregionalmax(amp);
  locMax1D = find(bw);
  [x,y,z]=ind2sub(size(amp),locMax1D);
  locs{i} = [x,y,z,amp(locMax1D)];
end

% Correct photobleach.
if license('test','Curve_Fitting_Toolbox')
  t = (0:nFrames-1)';
  pbFun = fit(t,meanInt,'exp1');
  pb0 = pbFun(t(1));
  for i=1:nFrames
    locs{i}(:,4) = pb0*locs{i}(:,4)/pbFun(t(i));
  end
else
  warning('Curve fitting toolbox unavailable. Not correcting for photobleach.');
end
if verbose
  figure;
  plot(t,meanInt,t,pbFun(t),t,pb0*meanInt./pbFun(t));
end

% Define bounds.
maxAmps = zeros(nFrames,1);
for i=1:nFrames
  maxAmps(i) = max(locs{i}(:,4));
end
minThresh = 0; % mean background level
maxThresh = min(maxAmps)-eps;

% Do global optimize.
opts = psoptimset('display','off','tolfun',1e-3,'cache','on','timelimit',300);
if verbose
  opts = psoptimset(opts,'outputfcns',@progress);
end
[threshold,fval] = patternsearch(@objective,repmat(0.5*(minThresh+maxThresh),1,nPieces), ...
    [],[],[],[],repmat(minThresh,1,nPieces),repmat(maxThresh,1,nPieces),[], opts);
% Go over all frames and apply the threshold.
spots = cell(nFrames,1);
for i=1:nFrames
  spots{i} = findSpots(locs{i}(:,1:3),locs{i}(:,4),piecewise_const_threshold(threshold,i,nPieces,nFrames));
if verbose
  figure;
  histogram(locs{i}(:,4));
  plot([piecewise_const_threshold(threshold,i,nPieces,nFrames) piecewise_const_threshold(threshold,i,nPieces,nFrames)],[0 1],'r--')
%   allAmps = cell2mat(locs);
%   allAmps = allAmps(:,4);
%   [f,x] = ksdensity(allAmps);
%   [~,xmax] = max(f);
%   plot(x,f,[threshold threshold],[0 max(f)],'r--');
end
end

% END


function [stop,options,optchanged] = progress(optimvalues,options,flag)
  stop = 0;
  optchanged = 0;
  n = 0;
  if ~strcmp(flag,'interrupt')
    for i=1:length(locs)
      n = n + sum(locs{i}(:,4)>=piecewise_const_threshold(optimvalues.x,i,nPieces,nFrames));
    end
    fprintf('Threshold: %.3f%%  Mean min diff: %g Mean spot count: %.1f\n',((optimvalues.x(1)+optimvalues.x(2)*i)-minThresh)/(maxThresh-minThresh),optimvalues.fval,n/length(locs));
  end
end

function thresh = piecewise_const_threshold(t,ind,nPieces,nFrames)
    %t: vector with length nPieces: takes nPieces different constant values throughout a movie
    %ind, nPieces, nFrames: all integers
    which_piece = idivide(int32(ind-1),ceil(nFrames/nPieces)) + 1; %work out which fraction of the movie we are in
    %which_piece should be between 1 and nPieces
    thresh = t(which_piece);
end

function y = objective(t)
    %objective to minimize as a function of the threshold value, t
  m = zeros(length(locs)-1,1);
  n = zeros(length(locs),1);
  for i=1:length(locs)-1
    % Find local maxima passing threshold in this frame and the next.
    s1 = findSpots(locs{i}(:,1:3),locs{i}(:,4),piecewise_const_threshold(t,i,nPieces,nFrames));
    s2 = findSpots(locs{i+1}(:,1:3),locs{i+1}(:,4),piecewise_const_threshold(t,i+1,nPieces,nFrames));
    % Compute metric for point cloud difference.
    m(i) = meanMinDiff(s1,s2);
    % Mean number of spots for optional penalty.
    n(i) = size(s1,1);
  end
  n(end) = size(s2,1);
  y = mean(m) + lambda*mean(abs(n-realisticNumSpots));
end

function s=findSpots(locMax,amp,th)
% Get local maxima from the image
  s = locMax(amp>=th,:);
end

end


