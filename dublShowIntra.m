function dublShowIntra(intraStruct,varargin)
% DUBLSHOWINTRA Plots both 3D and 2D scatters of outer-kinetochore position
% relative to the inner kinetochore marker.
%
%    DUBLSHOWINTRA(JOB,...) Plots a 3D scatter of outer-kinetochore
%    measurements in delta x, y and z coordinates relative to the average
%    inner-kinetochore marker as per the JOB. Also plots 2D scatters in
%    each the x-y and x-z planes. Rotates the data by default into the
%    plane of the inter-kinetochore axis, aligned along the x-axis. Can
%    also fit shells and circles with radii the median measurement of delta
%    3D. Other options available.
%
%    Options, defaults in {}:-
%
%    depthFilter: 0 or {1}. Whether or not to plot data filtered for delta
%           z between ±100nm.
%
%    fit: {0} or 1. Whether or not to fit a shell to the 3D data, which
%           will have radius the median of 3D delta.
%
%    rotate: 0 or {1}. Whether or not to rotate all data so that the x-axis
%           represents the inter-kinetochore axis. This will be represented
%           as a red line in all plots.
%
%    transparency: 0.5 or number in [0 1]. The transparency of the fitted
%           shell, so that 0 is entirely opaque, and 1 is invisible.
%
% Copyright (c) 2017 C. A. Smith

opts.depthFilter = 1;
opts.fit = 0;
opts.plotType = 'scatter'; % can also be 'heatMap'
opts.rotate = 1;
opts.transparency = 0.5;
opts = processOptions(opts,varargin{:});

% invert transparency
opts.transparency = 1-opts.transparency;

% get coordinate-specific delta
if opts.depthFilter
    delta(:,1) = intraStruct.plate.depthFilter.delta.y.all(:);
    delta(:,2) = intraStruct.plate.depthFilter.delta.x.all(:);
    delta(:,3) = intraStruct.plate.depthFilter.delta.z.all(:);
    delta3D    = intraStruct.plate.depthFilter.delta.threeD.all(:);
else
    delta(:,1) = intraStruct.plate.raw.delta.y.all(:);
    delta(:,2) = intraStruct.plate.raw.delta.x.all(:);
    delta(:,3) = intraStruct.plate.raw.delta.z.all(:);
    delta3D    = intraStruct.plate.raw.delta.threeD.all(:);
end

% and coordinate-specific sister separation
sisSep(:,1) = intraStruct.plate.sisSep.x;
sisSep(:,2) = intraStruct.plate.sisSep.y;
sisSep(:,3) = intraStruct.plate.sisSep.z;
sisSep = [sisSep;-sisSep];

% get x-positions
positX = intraStruct.plate.coords.x(:,[1 3]);
positX = [positX; fliplr(positX)];

% get numbers
nKTs = size(delta,1);

%% Process rotation of coordinates and fit

if opts.rotate
    
    % calculate angles for rotation
    % XY
    tanthe = sisSep(:,2)./sisSep(:,1);
    theta = atan(tanthe);
    % XZ
    tanrho = sisSep(:,3)./sisSep(:,1);
    rho = atan(tanrho);

    % loop over kinetochores to rotate them on K-K axis
    for iSis = 1:nKTs

        % get rotation angles for each x-y
        iTheta = theta(iSis);
        % and x-z
        iRho = rho(iSis);

        if positX(iSis,2)>positX(iSis,1)
          iTheta = iTheta+pi;
          iRho = iRho+pi;
        end
        % rotate delta
        delta(iSis,:) = delta(iSis,:)*[cos(iTheta),  sin(iTheta), 0;
                                       -sin(iTheta), cos(iTheta), 0;
                                       0,            0,           1];
        delta(iSis,:) = delta(iSis,:)*[cos(iRho),    0,           sin(iRho);
                                       0,            1,           0;
                                       -sin(iRho),   0,           cos(iRho)];
    end

end

%% Fitting a shell to the data

% remove NaNs
delta = delta(~isnan(delta(:,1)),:);

% calculate best fit shell
ft = fittype('r^2 - x^2 - y^2','coefficients',{'r'},'independent',{'x','y'},'dependent','z');
inputDel = nanmedian(delta3D); inputDel = [inputDel inputDel/10 inputDel*10];
shellFit = fit([delta(:,2),delta(:,1)],delta(:,3).^2,ft,...
    'StartPoint',inputDel(1),'Lower',inputDel(2),'Upper',inputDel(3),...
    'Robust','LAR');
fitR = coeffvalues(shellFit)*1000;
kitLog('Shell fitting yielded a measured radius of %.1f.',fitR);

% produce sphere for plotting
[sphX,sphY,sphZ] = sphere(20);
sphX = sphX*fitR; sphY = sphY*fitR; sphZ = sphZ*fitR;
% produce x-coordinates for circle plotting
th = 0:pi/50:2*pi;

%% Print figures

% ensure figures are docked
set(0,'DefaultFigureWindowStyle','docked')

% plot 3D scatter
figure;
scatter3(delta(:,2)*1000,delta(:,1)*1000,delta(:,3)*1000,'go',...
    'MarkerFaceColor','g')
% plot shell and inter-kinetochore axis
hold on
scatter3(0,0,0,'ro','MarkerFaceColor','r','LineWidth',3)
if opts.rotate; line([-150 0],[0 0],[0 0],'Color','r','LineWidth',2); end
if opts.fit
    surf(sphX,sphY,sphZ,...
        'FaceColor',[0 0.5 0],'FaceAlpha',opts.transparency,'EdgeColor','none')
end
% aesthetics
xlabel('\Delta_{x} (nm)'); ylabel('\Delta_{y} (nm)'); zlabel('\Delta_{z} (nm)');
xlim([-150 150]); ylim([-150 150]); zlim([-150 150]);
set(gca,'FontSize',20);

% plotting individual levels
incDelta = [-0.1 -0.07:0.02:0.07 0.1];

switch opts.plotType
    
    case 'scatter'
        
        figure;
        for i=1:9

            subplot(3,3,i);

            % produce scatter
            sel=(delta(:,3)>incDelta(i) & delta(:,3)<incDelta(i+1));
            scatter(delta(sel,2)*1000,delta(sel,1)*1000,'o',...
                'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);

            % plot shell and inter-kinetochore axis
            hold on
            scatter(0,0,'ro','MarkerFaceColor','r','LineWidth',3)
            if opts.rotate; line([-150 0],[0 0],[0 0],'Color','r','LineWidth',2); end
            if opts.fit && nanmean(incDelta(i:i+1)*1000) < fitR
                thisR = sqrt(fitR^2 - nanmean(incDelta(i:i+1)*1000)^2);
                xunit = thisR*cos(th); yunit = thisR*sin(th);
                plot(xunit,yunit,'Color',[0 0.5 0],'LineWidth',2)
            end

            % aesthetics
            xlabel('\Delta_{x} (nm)'); ylabel('\Delta_{y} (nm)');
            subTit = sprintf('%.0f < \\Delta_{z} < %.0f nm',incDelta(i)*1000,incDelta(i+1)*1000);
            title(subTit);
            xlim([-150 150]); ylim([-150 150]);

        end

        figure;
        for i=1:9

            subplot(3,3,i);

            % produce scatter
            sel=(delta(:,1)>incDelta(i) & delta(:,1)<incDelta(i+1));
            scatter(delta(sel,2)*1000,delta(sel,3)*1000,'o',...
                'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);

            % plot shell and inter-kinetochore axis
            hold on
            scatter(0,0,'ro','MarkerFaceColor','r','LineWidth',3)
            if opts.rotate; line([-150 0],[0 0],[0 0],'Color','r','LineWidth',2); end
            if opts.fit && nanmean(incDelta(i:i+1)*1000) < fitR
                thisR = sqrt(fitR^2 - nanmean(incDelta(i:i+1)*1000)^2);
                xunit = thisR*cos(th); yunit = thisR*sin(th);
                plot(xunit,yunit,'Color',[0 0.5 0],'LineWidth',2)
            end

            % aesthetics
            xlabel('\Delta_{x} (nm)'); ylabel('\Delta_{z} (nm)');
            subTit = sprintf('%.0f < \\Delta_{y} < %.0f nm',incDelta(i)*1000,incDelta(i+1)*1000);
            title(subTit);
            xlim([-150 150]); ylim([-150 150]);

        end
        
    case 'heatMap'
        
        figure;
        for i=1:9

            subplot(3,3,i);

            % produce scatter
            sel=(delta(:,3)>incDelta(i) & delta(:,3)<incDelta(i+1));
            heatMap(delta(sel,2)*1000,delta(sel,1)*1000,'nBins',15,...
                'xLimits',[-150 150],'yLimits',[-150 150],'withinFig',1)

            % plot shell and inter-kinetochore axis
            hold on
            scatter(0,0,'wx','MarkerFaceColor','w','LineWidth',2)
            if opts.rotate; line([-150 0],[0 0],[0 0],'Color','w','LineWidth',2); end
            if opts.fit && nanmean(incDelta(i:i+1)*1000) < fitR
                thisR = sqrt(fitR^2 - nanmean(incDelta(i:i+1)*1000)^2);
                xunit = thisR*cos(th); yunit = thisR*sin(th);
                plot(xunit,yunit,'Color','w','LineWidth',2)
            end

            % aesthetics
            xlabel('\Delta_{x} (nm)'); ylabel('\Delta_{y} (nm)');
            subTit = sprintf('%.0f < \\Delta_{z} < %.0f nm',incDelta(i)*1000,incDelta(i+1)*1000);
            title(subTit);

        end

        figure;
        for i=1:9

            subplot(3,3,i);

            % produce scatter
            sel=(delta(:,1)>incDelta(i) & delta(:,1)<incDelta(i+1));
            heatMap(delta(sel,2)*1000,delta(sel,3)*1000,'nBins',15,...
                'xLimits',[-150 150],'yLimits',[-150 150],'withinFig',1)

            % plot shell and inter-kinetochore axis
            hold on
            scatter(0,0,'wx','MarkerFaceColor','w','LineWidth',2)
            if opts.rotate; line([-150 0],[0 0],[0 0],'Color','w','LineWidth',2); end
            if opts.fit && nanmean(incDelta(i:i+1)*1000) < fitR
                thisR = sqrt(fitR^2 - nanmean(incDelta(i:i+1)*1000)^2);
                xunit = thisR*cos(th); yunit = thisR*sin(th);
                plot(xunit,yunit,'Color','w','LineWidth',2)
            end

            % aesthetics
            xlabel('\Delta_{x} (nm)'); ylabel('\Delta_{z} (nm)');
            subTit = sprintf('%.0f < \\Delta_{y} < %.0f nm',incDelta(i)*1000,incDelta(i+1)*1000);
            title(subTit);

        end
    otherwise
        
        warning('Plot type not recognised: %s.',opts.plotType);
        
end

end



