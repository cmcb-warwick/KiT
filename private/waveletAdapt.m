function tk = waveletAdapt(movie,options)
% Compare successive frames on nSpots and covaraince of distance matrix difference

kitLog('Determining adaptive wavelet threshold');

% tk = 0.1;
% tkinc = 0.1;
% tkincfac = 1.1;
% tkmax = 50;
% %f = [1 floor(nFrames/4) floor(nFrames/2) floor(3*nFrames/4) nFrames-1];
nFrames = size(movie,4);
% f = [1  floor(nFrames/2) nFrames-1];
% %f = [1 nFrames-1];
% first = 1;
% i = 1;
% while first || (tk < tkmax && nSpots > 0 && ~isnan(ld))
%   options.waveletLevelThresh=tk;
%   for k=1:length(f)
%     [A,~,ld] = waveletSpots(movie(:,:,:,f(k)),options);
%     B = waveletSpots(movie(:,:,:,f(k)+1),options);
%     frameDiff(i,k) = meanMinDiff(A,B);
%   end

%   % Increment tk.
%   first = 0;
%   tkvec(i) = tk;
%   tk = tk + tkinc;
%   tkinc = tkinc*tkincfac;
%   nSpots = size(A,1);
%   i = i+1;
% end
% % Pick tk which minimises frameDiff metric, with small penalty for increasing tk.
% lambda = 0.01; % penalty factor.
% frameDiff = mean(frameDiff,2) + lambda*tkvec';
% [~,minIdx] = min(frameDiff);
% pp = pchip(tkvec,frameDiff);
% tk = fminbnd(@(x) ppval(pp,x),tkvec(max(minIdx-1,1)),tkvec(min(minIdx+1,length(tkvec))));
lambda = 10;
opts = psoptimset('display','off','tolfun',1e-3,'cache','on','timelimit',300,'plotfcns',{@psplotbestf,@psplotbestx});
tk = patternsearch(@objective,2,[],[],[],[],1,50,[], opts);


kitLog('Using wavelet threshold: %g',tk);

if options.debug.showWaveletAdapt
  % figure;
  % x = linspace(tkvec(1),tkvec(end));
  % plot(x,ppval(pp,x));
  % hold on
  % plot([tk tk],ylim,'r--');
  % hold off
  % ylabel('Frame-frame point cloud difference');
  % xlabel('Wavelet threshold');
  % drawnow;
end


function y = objective(t)
  f = [1  floor(nFrames/2) nFrames-1];
  frameDiff = zeros(length(f),1);
  options.waveletLevelThresh = t;
  for k=1:length(f)
    A = waveletSpots(movie(:,:,:,f(k)),options);
    B = waveletSpots(movie(:,:,:,f(k)+1),options);
    if isempty(A) || isempty(B)
      m(k,1) = inf;
    else
      % Compute metric for point cloud difference.
      m(k,1) = meanMinDiff(A,B);
    end
    % Mean number of spots for optional penalty.
    m(k,2) = 0.5*size(A,1)*size(B,1);
  end
  %y = mean(m(:,1)) + lambda*t;
  y = mean(m(:,1)) + lambda/(max(1,mean(m(:,2))^.25));
end


end
