function planeFit = smoothOriginTimeSeries(planeFit, verbose)
%% Smooth the jumpy movement of the coordinate origin.
%% These movements may be partially due to
%%%%%%%%%%%%%%%%%%%%
nFrames = length(planeFit);
probDim = size(planeFit(1).planeOrigin,2); 
originCoords = zeros(nFrames,3);
for j=1:nFrames
    originCoords(j,:) = planeFit(j).planeOrigin;
end
smoothedOriginCoords = smoothdata(originCoords,'movmean',10);
for j=1:nFrames
    planeFit(j).planeOrigin = smoothedOriginCoords(j,:);
end
if verbose
    figure;
    for i=1:probDim        
        subplot(probDim,1,i);
        plot(originCoords(:,i),'linewidth',3);
        hold all;
        plot(smoothedOriginCoords(:,i),'linewidth',3);
        xlabel('Frame');
        ylabel(sprintf('Position, um'));
        set(gca,'fontsize',20);
    end
end
