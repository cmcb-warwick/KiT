function [spots,amps,bgAmps,rejects]=mixtureModelFit(cands,image,psfSigma,options)
% MIXTUREMODELFIT Fit Gaussian mixture models to spots
%
% Enhances candidate spots obtained by centroid fitting by fitting Gaussian
% mixture models. F-test is used to determine appropriate number of
% Gaussians.
%
% Copyright (c) 2014 J. W. Armond

verbose = options.debug.mmfVerbose;

% Set optimization options.
optoptions = optimset('Jacobian','on','Display','off','Tolfun',1e-4,'TolX',1e-4);

alphaF = options.alphaF; % N vs N+1 F-test cutoff.
alphaA = options.alphaA; % amplitude t-test cutoff.
alphaD = options.alphaD;
is3D = ndims(image) == 3;
if is3D
  degFreePerSpot = 4;
  cols = 3;
  psfSigma = psfSigma([1 1 2]); % for XYZ
else
  degFreePerSpot = 3;
  cols = 2;
  psfSigma = psfSigma([1 1]); % for XY
end
% Distance criterion for defining separate clusters
clusterSep = options.clusterSeparation*psfSigma;

% Cluster spots.
clusters = clusterSpots(cands(:,1:cols));
nClusters = max(clusters);
if verbose
  kitLog('Mmf fitting on %d clusters',nClusters);
end
spots = [];
amps = [];
bgAmps = [];
rejects.amp = []; % Candidates rejected on amplitude test.
rejects.dist = []; % Candidate rejected on distance test.

% Disable singular matrix warning. Some clusters are ill-conditioned and cause
% the warning in the variance, but is harmless
warning('off','MATLAB:singularMatrix');

startTime = clock;

  % For each cluster, perform iterative mixture-model fitting, increasing
  % number of Gaussians and F-testing.
  for i=1:nClusters
    % Extract cluster data.
    clusterCandsT = cands(clusters==i,:); % Candidates in this cluster.
    if is3D
      clusterCands1D = sub2ind(size(image),clusterCandsT(:,1),clusterCandsT(:,2),...
                               clusterCandsT(:,3));
    else
      clusterCands1D = sub2ind(size(image),clusterCandsT(:,1),clusterCandsT(:,2));
    end
    clusterAmpT = image(clusterCands1D);
    numCandsT = size(clusterCandsT,1);
    if verbose
      kitLog('Fitting cluster %d, %d cands',i,numCandsT);
    end

    clusterPix = getClusterPixels(clusterCandsT(:,1:cols)); % Pixels in this cluster.
    if is3D
      clusterPix1D = sub2ind(size(image),clusterPix(:,1),clusterPix(:,2),...
                             clusterPix(:,3));
    else
      clusterPix1D = sub2ind(size(image),clusterPix(:,1),clusterPix(:,2));
    end
    clusterImg = image(clusterPix1D);

    % Background amplitude estimate, initial guess.
    bgAmpT = mean(clusterImg(:));

    % Iterative fitting with N+1 Gaussians until F-test fails.
    first = 1; % Always accept first fit.
    failed = 0;
    while (~failed && options.mmfAddSpots) || first
      % Estimate bounds.
      [x0,lb,ub] = guessBounds(clusterCandsT(:,1:cols),clusterAmpT,clusterPix,bgAmpT);

      % Fit mixture-model.
      [solutionT,resnorm,residuals,jacT] = fitNGaussiansFitFun(optoptions,x0,lb,ub,...
                                                        clusterImg,clusterPix,psfSigma);

      numDegFreeT = size(clusterPix,1) - degFreePerSpot*numCandsT - 1;
      residVarT = resnorm/numDegFreeT;

      if first
        if verbose
          kitLog('First fit complete');
        end
        first = 0;
      else
        % F distributed test statistic.
        testStat = residVarT/residVar;
        pValue = fcdf(testStat,numDegFree,numDegFreeT);
        if pValue > alphaF
          failed = 1;
        end
        if verbose
          if failed
            resultStr = 'rejected';
          else
            resultStr = 'accepted';
          end
          kitLog('N+1 fit complete, %d cands, p=%f (%s)',numCandsT,pValue,resultStr);
        end
      end % if first

      if ~failed
        % Update accepted variables.
        numCands = numCandsT;
        numDegFree = numDegFreeT;
        solution = solutionT;
        residVar = residVarT;
        jac = jacT;

        % Extract values from solution vector.
        [bgAmp,clusterCands,clusterAmp] = extractSolution(solution,numCands);

        if options.mmfAddSpots
          % Add new kernel at pixel with maximum residual.
          numCandsT = numCandsT + 1;
          clusterAmpT = [clusterAmp; mean(clusterAmp);];
          [~,idx] = max(residuals);
          coord = clusterPix(idx,:);
          clusterCandsT = [clusterCands; coord];
          bgAmpT = bgAmp;
        end
      end % if pValue > alphaF

    end % while ~failed

    if options.maxMmfTime > 0 && etime(clock,startTime) > options.maxMmfTime
      spots = []; amps = []; bgAmps = [];
      return % Abort
    end

    % Accumulate spots.
    spots = [spots; clusterCands];
    amps = [amps; clusterAmp];
    bgAmps = [bgAmps; repmat(bgAmp,[numCands,1])];
  end

% Recluster to incorporate new candidates that are nearby into same
% cluster.
if isempty(spots)
  return
end
clusters = clusterSpots(spots(:,1:cols));
nClusters = max(clusters);
if verbose
    kitLog('Reclustering produced %d clusters',nClusters);
end

% New output variables.
spots2 = [];
amps2 = [];
bgAmps2 = [];

% Distance testing. 1-sided t-test.
for i=1:nClusters
  % Extract candidate information for new cluster.
  idx = clusters==i;
  clusterCands = spots(idx,1:cols);
  clusterAmp = amps(idx,1);
  firstIdx = find(idx,1);
  bgAmp = bgAmps(firstIdx,1);
  numCands = size(clusterCands,1);
  assert(numCands > 0);
  clusterPix = getClusterPixels(clusterCands); % Pixels in this cluster.
  if is3D
    clusterPix1D = sub2ind(size(image),clusterPix(:,1),clusterPix(:,2),...
                           clusterPix(:,3));
  else
    clusterPix1D = sub2ind(size(image),clusterPix(:,1),clusterPix(:,2));
  end
  clusterImg = image(clusterPix1D);

  % Refit.
  % Estimate bounds.
  [x0,lb,ub] = guessBounds(clusterCands(:,1:cols),clusterAmp,clusterPix,bgAmp);

  % Refit mixture-model.
  [solution,resnorm,~,jac] = fitNGaussiansFitFun(optoptions,x0,lb,ub,...
                              clusterImg,clusterPix,psfSigma);

  numDegFree = size(clusterPix,1) - degFreePerSpot*numCands - 1;
  residVar = resnorm/numDegFree;

  % Extract values from solution vector.
  [bgAmp,clusterCands,clusterAmp] = extractSolution(solution,numCands);
  if verbose
      kitLog('Refit to new cluster, %d cands',numCands);
  end

  while numCands > 1
    % Calculate p-values of distances between all candidates.
    [~,~,~,covMat] = computeVariances(jac,residVar,numCands);
    pValue = mmfDistPV(clusterCands,covMat,numCands,numDegFree);

    % Maximum p-value.
    [pValueMax,indxBad] = max(pValue(:));
    testDist = pValueMax > alphaD;
    if ~testDist
      if verbose
        kitLog('All candidates passed distance test. Range of p=[%g,%g]', ...
               min(pValue(:)),pValueMax);
      end
      break;
    end

    % Identify pair with maximum p-value.
    [indx1,indx2] = ind2sub(size(pValue),indxBad);
    % Identify candidate with smaller amplitude from pair.
    ampPair = clusterAmp([indx1 indx2]);
    if ampPair(1) < ampPair(2)
      indxBad = indx1;
    else
      indxBad = indx2;
    end

    % Remove candidate.
    rejects.dist = [rejects.dist; clusterCands(indxBad,:) pValueMax];
    clusterCands(indxBad,:) = [];
    clusterAmp(indxBad,:) = [];
    numCands = numCands - 1;
    if numCands == 0
        if verbose
          kitLog('All candidates removed due to failing distance test');
        end
        break;
    end

    % Estimate bounds.
    [x0,lb,ub] = guessBounds(clusterCands(:,1:cols),clusterAmp,clusterPix,bgAmp);

    % Refit mixture-model.
    [solution,resnorm,~,jac] = fitNGaussiansFitFun(optoptions,x0,lb,ub,...
                                   clusterImg,clusterPix,psfSigma);

    numDegFree = size(clusterPix,1) - degFreePerSpot*numCands - 1;
    residVar = resnorm/numDegFree;

    % Extract values from solution vector.
    [bgAmp,clusterCands,clusterAmp] = extractSolution(solution,numCands);

    if verbose
      kitLog('Removed candidate failing distance test, p=%f, %d cands',pValueMax,numCands);
    end

  end

  if numCands == 0
    continue;
  end

  % Compute variance estimates.
  [~,clusterAmpVar] = computeVariances(jac,residVar,numCands);

  % Amplitude testing. 1-sided t-test.
  % Repeat while some amplitudes not signifcant.
  while numCands > 0
    testStat = clusterAmp./sqrt(clusterAmpVar+residVar);
    pValue = 1-tcdf(testStat,numDegFree);
    [pValueMax,indxBad] = max(pValue);
    testAmp = pValueMax > alphaA;
    if ~testAmp
      if verbose
        kitLog('All candidates passed amplitude test. Range of p=[%g,%g]',min(pValue),pValueMax);
      end
      break;
    end

    % Remove candidate.
    rejects.amp = [rejects.amp; clusterCands(indxBad,:) clusterAmp(indxBad,:) ...
                  pValueMax];
    clusterCands(indxBad,:) = [];
    clusterAmp(indxBad,:) = [];
    numCands = numCands - 1;
    if numCands == 0
      if verbose
        kitLog('All candidates removed due to failing ampltitude test');
      end
      break;
    end

    % Estimate bounds.
    [x0,lb,ub] = guessBounds(clusterCands(:,1:cols),clusterAmp,clusterPix,bgAmp);

    % Refit mixture-model.
    [solution,resnorm,~,jac] = fitNGaussiansFitFun(optoptions,x0,lb,ub,...
                                   clusterImg,clusterPix,psfSigma);

    numDegFree = size(clusterPix,1) - degFreePerSpot*numCands - 1;
    residVar = resnorm/numDegFreeT;

    % Extract values from solution vector.
    [bgAmp,clusterCands,clusterAmp] = extractSolution(solution,numCands);
    % Compute variance estimates.
    [~,clusterAmpVar] = computeVariances(jac,residVar,numCands);

    if verbose
      kitLog('Removed candidate failing ampltitude test, p=%f, %d cands',pValueMax,numCands);
    end
  end

  if numCands == 0
    continue;
  end

  % Compute variance estimates.
  [clusterCandsVar,clusterAmpVar,bgAmpVar] = computeVariances(jac,residVar,numCands);

  if verbose
    kitLog('Cluster fitting complete, %d spots',numCands);
  end

  % Add amplitude p-value
  testStat = clusterAmp./sqrt(clusterAmpVar+residVar);
  pValue = 1-tcdf(testStat,numDegFree);

  if options.maxMmfTime > 0 && etime(clock,startTime) > options.maxMmfTime
    spots = []; amps = []; bgAmps = [];
    return % Abort
  end

  % Convert clusters into array of spots.
  spots2 = [spots2; [clusterCands clusterCandsVar]];
  amps2 = [amps2; [clusterAmp clusterAmpVar pValue]];
  bgAmps2 = [bgAmps2; repmat([bgAmp bgAmpVar],[numCands,1])];
end

spots = spots2;
amps = amps2;
bgAmps = bgAmps2;

if ~isempty(spots)
  % Calculate standard deviations from variances.
  spots(:,cols+1:end) = sqrt(spots(:,cols+1:end));
  amps(:,2) = sqrt(amps(:,2));
  bgAmps(:,2) = sqrt(bgAmps(:,2));

  % Convert to image coordinates.
  if is3D
    spots = spots(:,[2 1 3 5 4 6]);
  else
    spots = spots(:,[2 1 4 3]);
  end
end

if ~isempty(rejects.amp)
  rejects.amp(:,1:2) = rejects.amp(:,[2 1]);
end
if ~isempty(rejects.dist)
  rejects.dist(:,1:2) = rejects.dist(:,[2 1]);
end

% Re-enable warnings.
warning('on','MATLAB:singularMatrix');

if verbose
  kitLog('MMF complete, %d spots',size(spots,1));
end


%% SUBFUNCTIONS

function c=clusterSpots(pos)
% pos: spot positions (nx3)

  if options.oneBigCluster~=0
    c = ones(size(pos,1),1);
    return;
  end

  % Only 1 spot?
  if size(pos,1)==1
    c = 1;
  else
    % Form hierarchical clusters. Standardize Euclidean distances by dividing
    % by separation factor.
    z = linkage(pos,'single',{'seuclidean',clusterSep});
    c = cluster(z,'cutoff',1,'criterion','distance');
    c = c(:,1);
  end

  % Visualization.
  if options.debug.showMmfClusters ~= 0

    % If 3D image, max project.
    img = max(image,[],3);

    % Make 3 layers out of original image (normalized)
    img = img/max(img(:));
    img = repmat(img,[1 1 3]);

    % Get set of colors for labeling
    nClusters = max(c);
    color = jet(nClusters);
    color = color(randperm(nClusters),:);

    % Label maxima (each island has a different color)
    pix = round(pos(:,1:2));
    for j=1:size(pos,1)
        img(pix(j,1),pix(j,2),:)=color(c(j),:);
    end

    %plot image
    figure(1);
    imshow(img);
    title('Overlapping PSF clusters');
    drawnow;
    switch options.debug.showMmfClusters
      case -1
        pause;
      case -2
        keyboard;
    end

    % scatter version
    % scatter3(pos(:,1),pos(:,2),pos(:,3),50,c(:,1),'filled')
  end

end % function clusterSpots


function [x0,lb,ub]=guessBounds(pos,amp,clusterPixels,bgAmp)
% pos: candidate maxima positions (nx3)
% clusterPixels: coordinates of pixels comprising the cluster (nx3)

  x0 = pos; % initial guess is candidate positions

  % Lower bound.
  lb = x0;
  lb(:,1:2) = lb(:,1:2) - 2*max(1,psfSigma(1));
  if is3D
    lb(:,3) = lb(:,3) - max(1,psfSigma(3));
  end
  minPos = min(clusterPixels,[],1); % Bound by pixels comprising cluster
  lb = bsxfun(@max,lb,minPos);

  % Upper bound.
  ub = x0;
  ub(:,1:2) = ub(:,1:2) + 2*max(1,psfSigma(1));
  if is3D
    ub(:,3) = ub(:,3) + max(1,psfSigma(3));
  end
  maxPos = max(clusterPixels,[],1); % Bound by pixels comprising cluster
  ub = bsxfun(@min,ub,maxPos);

  % Add amplitudes.
  x0 = [x0 amp];
  lb(:,cols+1) = eps;
  ub(:,cols+1) = 1;

  % Reshape and add background.
  x0 = x0';
  x0 = [x0(:); bgAmp];
  lb = lb';
  lb = [lb(:); eps];
  ub = ub';
  ub = [ub(:); 1];

end % function guessBounds


function pixelCoords=getClusterPixels(clusterPos,visual)
% clusterPos: coordinates of pixels in cluster (nx3)
  if nargin<2
    visual = 0;
  end

  bound = clusterSep;
  minPos = max(1,floor(min(clusterPos,[],1)-bound));
  maxPos = min(size(image),ceil(max(clusterPos,[],1)+bound));

  if is3D
    [x,y,z] = ndgrid(minPos(1):maxPos(1),minPos(2):maxPos(2),...
                     minPos(3):maxPos(3));
    pixelCoords = [x(:) y(:) z(:)];
  else
    [x,y] = ndgrid(minPos(1):maxPos(1),minPos(2):maxPos(2));
    pixelCoords = [x(:) y(:)];
  end

  if visual
    figure(1);
    % If 3D image, max project.
    img = max(image(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)),[],3);

    % Make 3 layers out of original image (normalized)
    img = img/max(img(:));
    img = repmat(img,[1 1 3]);

    % Label cluster spots
    pix = round(bsxfun(@minus,clusterPos(:,1:2),minPos(1:2))+1);
    for j=1:size(pix,1)
        img(pix(j,1),pix(j,2),:)=[1 0 0];
    end
    imshow(img);
  end
end % function pixelCoords

function [bgAmp,clusterCands,clusterAmp]=extractSolution(solution,numCands)
% solution: output from fitting.

  % Background amplitude estimate.
  bgAmp = solution(end);
  solution(end) = [];

  % Reshape solution in nx4.
  solution = reshape(solution,cols+1,numCands)';

  % Extract positions and amplitudes.
  clusterCands = solution(:,1:cols);
  clusterAmp = solution(:,cols+1);
end

function [candVar,ampVar,bgAmpVar,covMat]=computeVariances(jac,residVar,numCands)
  jac = full(jac);
  jac = jac'*jac;

  % Estimate covariances from Jacobian.
  if rcond(jac) < 1e-12
    covMat = pinv(jac)*residVar;
  else
    covMat = inv(jac)*residVar;
  end
  sigmaSq = diag(covMat);

  % Background amplitude variance.
  bgAmpVar = sigmaSq(end);
  sigmaSq(end) = [];

  % Reshape variance vector to nx4.
  sigmaSq = reshape(sigmaSq,cols+1,numCands)';

  % Extract position and amplitude variances.
  candVar = sigmaSq(:,1:cols);
  ampVar = sigmaSq(:,cols+1);
end

function showSpots(spots,z)
  if nargin<2
    z = 0;
  end

  figure(1);
  if z == 0
    % If 3D image, max project.
    img = max(image,[],3);
  else
    img = image(:,:,z);
  end

  % Make 3 layers out of original image (normalized)
  img = img/max(img(:));
  img = repmat(img,[1 1 3]);

  % Label cluster spots
  pix = round(spots);
  for j=1:size(pix,1)
    if z == 0 || pix(j,3) == z
      img(pix(j,1),pix(j,2),:)=[1 0 0];
    end
  end
  imshow(img);
end

end % function mixtureModelFit
