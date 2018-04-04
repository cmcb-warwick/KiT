function analysis = chrsAnalysis(varargin)
% CHRSANALYSIS Overall analysis of chromatic shift for multiple experiments.
%
%    CHRSANALYSIS(...) Produces figures of multiple chromatic shift
%    measurements for analysis. Either the images or a pre-processed
%    zetaMeasurements structure can be provided. Options are available.
%
%    Options, defaults in {}:-
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. Only applies when
%       movieStructure provided. The channels between which to make
%       chromatic shift measurements. The direction of measurements will
%       be defined from the first to the second stated channel number,
%       e.g. [2 1] means from channel 2 to 1.
%
%    filename: {[]} or string. Name of experiment, which will be present
%       in all output files.
%
%    zetaStructure: {[]} or output from chrsZetaMeasurements. Analyses
%       will be drawn from the zetaMeasurements structure output from
%       chrsZetaMeasurements.
%
%    imageStructure: {[]} or a cell of experiments. Analyses will be drawn
%       directly from the images across all experiments.
%
%    save: {0} or 1. Whether or not to save all figures to file. The user
%       will be asked to provide a save directory.
%
%    savePath: {[]} or string of directory path. Location in which all
%       analysis files should be saved.
%
%
% Copyright (c) 2017 C. A. Smith

% default options
opts.channels = [1 2];
opts.filename = [];
opts.zetaStructure = [];
opts.imageStructure = [];
opts.save = 0;
opts.savePath = [];
% process user-defined options
opts = processOptions(opts,varargin{:});

kitLog('Running chrsAnalysis')

zM = opts.zetaStructure;
iS = opts.imageStructure;

% if no intra- or movie structures provided, ask user to find one
if isempty(zM) && isempty(iS)
  % ask here for one or the other
  error('No zeta-measurements or image structure provided.')
end
% if movies provided, get intraMeasurements structure
if ~isempty(iS)
  zM = chrsZetaMeasurements(iS,'channels',opts.channels);
  opts.type = 'image structure';
else
  opts.type = 'zeta-measurements';
end

% check all requirements for saving
if opts.save
  
  % get a filename if not already provided
  while isempty(opts.filename)
    opts.filename = input('Please provide a file name: ','s');
  end
  
  % get a save directory if not already provided
  while isempty(opts.savePath)
    [~,opts.savePath] = uiputfile('*.mat','Save directory for all analysis files',[opts.filename '.mat']);
  end
    
end 
    
% set up default analysis structure
analysis = defaultAnalStruct;
analysis.filename = opts.filename;
analysis.filepath = opts.savePath;

%% Zeta XYZ

kitLog('Processing figures of coordinate-specific zeta');

% get measurements
zetaXYZ(:,1) = zM.filtered.zeta.x(:);
zetaXYZ(:,2) = zM.filtered.zeta.y(:);
zetaXYZ(:,3) = zM.filtered.zeta.z(:);

% find stats
mean_zetaXYZ = nanmean(zetaXYZ);
median_zetaXYZ = nanmedian(zetaXYZ);
std_zetaXYZ = nanstd(zetaXYZ);
serr_zetaXYZ = nanserr(zetaXYZ);
% store statistical information
analysis.stats.zetaX.mean   = mean_zetaXYZ(1);
analysis.stats.zetaY.mean   = mean_zetaXYZ(2);
analysis.stats.zetaZ.mean   = mean_zetaXYZ(3);
analysis.stats.zetaX.median = median_zetaXYZ(1);
analysis.stats.zetaY.median = median_zetaXYZ(2);
analysis.stats.zetaZ.median = median_zetaXYZ(3);
analysis.stats.zetaX.stdDev = std_zetaXYZ(1);
analysis.stats.zetaY.stdDev = std_zetaXYZ(2);
analysis.stats.zetaZ.stdDev = std_zetaXYZ(3);
analysis.stats.zetaX.stdErr = serr_zetaXYZ(1);
analysis.stats.zetaY.stdErr = serr_zetaXYZ(2);
analysis.stats.zetaZ.stdErr = serr_zetaXYZ(3);

% produce plot
chrsBasicPlots(zM,'stat','zetaXYZ','filtered',1,'legend','off')
if opts.save
  
  f = gcf;
  
  % create filename
  filename = ['hist_zetaXYZ_' opts.filename];
  % print the figure to file
  print(fullfile(opts.savePath,filename),'-depsc');
  
  close(f);
  
end

%% Positional correlations

kitLog('Processing figures of position-dependent zeta');

% get coordinates
coordsXYZ(:,1) = zM.coords.x(:);
coordsXYZ(:,2) = zM.coords.y(:);
coordsXYZ(:,3) = zM.coords.z(:);

% find the data for which there is information
sel = ~isnan(zetaXYZ(:,1));
analysis.stats.n = sum(sel);

% find all correlation coefficients
R = zeros(3);
p = ones(3);
coordLabels = 'xyz';
capsLabels = 'XYZ';
for i=1:3
  for j=1:3
    
    % find correlation
    [R(i,j),p(i,j)] = corr(coordsXYZ(sel,i),zetaXYZ(sel,j),'rows','pairwise','type','Spearman');
    
    % construct figure
    title = 'Position-dependent chromatic shift';
    xLab = ['position in ' coordLabels(i) ' (\mum)'];
    yLab = ['\zeta_{' coordLabels(j) '} (nm)'];
    xRange = [0 max(coordsXYZ(:,i))];
    yRange = [median_zetaXYZ(:,j)-1.5*std_zetaXYZ(:,j) median_zetaXYZ(:,j)+1.5*std_zetaXYZ(:,j)]*1000;
    overlay{1} = sprintf('\\rho = %.2f',R(i,j));
    overlay{2} = sprintf('p = %.3f',p(i,j));
    heatMap(coordsXYZ(sel,i),zetaXYZ(sel,j)*1000,'nBins',15,'withinFig',0,...
        'xLimits',xRange,'yLimits',yRange,'xLabel',xLab,'yLabel',yLab,'title',title);
    text(xRange(2)-diff(xRange)/20,yRange(2)-diff(yRange)/10,overlay,...
        'Color','w','FontSize',15,'HorizontalAlignment','right');
    
    if opts.save

      f = gcf;

      % create filename
      filename = ['heatMap_zeta' capsLabels(j) 'vs' capsLabels(i) '_' opts.filename];
      % print the figure to file
      print(fullfile(opts.savePath,filename),'-depsc');
  
      close(f);

    end
  end
end
% store information
cortn.R = R; cortn.p = p;
analysis.correlation = cortn;


%% Line of best fit

kitLog('Processing lines of best fit for position-dependent zeta');

[p1,p2] = find(p<0.01); p = [p1 p2]';
lineFit = analysis.lineFit;
if size(p,2)>0
  for i=p
    
    % fit polynomial-1 here for each i(1) and i(2)
    [coeffs,goodnss] = fit(coordsXYZ(sel,i(1)),zetaXYZ(sel,i(2))*1000,'poly1','Robust','Bisquare');
    coeffs = coeffvalues(coeffs);
    lineFit.ax(i(1),i(2)) = coeffs(1);
    lineFit.b(i(1),i(2)) = coeffs(2);
    lineFit.R2(i(1),i(2)) = goodnss.rsquare;
    
    % construct figure
    title = 'Linear approximation of position-dependent chromatic shift';
    xLab = ['position in ' coordLabels(i(1)) ' (\mum)'];
    yLab = ['\zeta_{' coordLabels(i(2)) '} (nm)'];
    overlay{1} = sprintf('\\zeta_{%s} = %.2f%s + %.2f',coordLabels(i(2)),coeffs(1),coordLabels(i(1)),coeffs(2));
    overlay{2} = sprintf('R^{2} = %.2f',goodnss.rsquare);
    xRange = [0 max(coordsXYZ(:,i(1)))];
    yRange = [median_zetaXYZ(:,i(2))-1.5*std_zetaXYZ(:,i(2)) median_zetaXYZ(:,i(2))+1.5*std_zetaXYZ(:,i(2))]*1000;
    heatMap(coordsXYZ(sel,i(1)),zetaXYZ(sel,i(2))*1000,'nBins',15,'withinFig',0,...
        'xLimits',xRange,'yLimits',yRange,...
        'xLabel',xLab,'yLabel',yLab,'title',title);
    hold on
    line(xRange,coeffs(1)*xRange+coeffs(2),'Color','w','LineWidth',2);
    text(xRange(2)-diff(xRange)/20,yRange(2)-diff(yRange)/10,overlay,...
        'Color','w','FontSize',15,'HorizontalAlignment','right');
    
    if opts.save

      f = gcf;

      % create filename
      filename = ['bestFitLine_zeta' capsLabels(i(2)) 'vs' capsLabels(i(1)) '_' opts.filename];
      % print the figure to file
      print(fullfile(opts.savePath,filename),'-depsc');
  
      close(f);

    end
    
  end
end
analysis.lineFit = lineFit;

if opts.save
  % save the results
  chrsAnalysis = analysis;
  save(fullfile(opts.savePath,opts.filename),'chrsAnalysis');
  kitLog('%s.mat and all figures saved in %s.',analysis.filename,analysis.filepath(1:end-1));
end

% print the results
printStatistics(analysis,opts);
kitLog('Analysis complete.');


end

%% Sub-functions

function printStatistics(aS,opts)

% output general statistics
kitLog('Summary of analysis for %s: %s...',opts.type,opts.filename);

stats = aS.stats;
kitLog('Distance measurements for zeta');
fprintf('\nPositional zeta [x,y,z] (n = %i):\n',stats.n)
fprintf('Median:  [%.1f, %.1f, %.1f] nm\n',...
    stats.zetaX.median*1000, stats.zetaY.median*1000, stats.zetaZ.median*1000);
fprintf('Mean:    [%.1f, %.1f, %.1f] nm\n',...
    stats.zetaX.mean*1000, stats.zetaY.mean*1000, stats.zetaZ.mean*1000);
fprintf('Std dev: [%.1f, %.1f, %.1f] nm\n',...
    stats.zetaX.stdDev*1000, stats.zetaY.stdDev*1000, stats.zetaZ.stdDev*1000);
fprintf('Std err: [%.1f, %.1f, %.1f] nm\n',...
    stats.zetaX.stdErr*1000, stats.zetaY.stdErr*1000, stats.zetaZ.stdErr*1000);
fprintf('\n')

% specific correlation information
cortn = aS.correlation;
kitLog('Results of correlation analysis');
fprintf('\nCorrelation coefficients:\n');
fprintf('Zeta x vs. x position:   rho = %.3f, p = %.3f\n', cortn.R(1,1), cortn.p(1,1));
fprintf('           y position:   rho = %.3f, p = %.3f\n', cortn.R(2,1), cortn.p(2,1));
fprintf('           z position:   rho = %.3f, p = %.3f\n', cortn.R(3,1), cortn.p(3,1));
fprintf('Zeta y vs. x position:   rho = %.3f, p = %.3f\n', cortn.R(1,2), cortn.p(1,2));
fprintf('           y position:   rho = %.3f, p = %.3f\n', cortn.R(2,2), cortn.p(2,2));
fprintf('           z position:   rho = %.3f, p = %.3f\n', cortn.R(3,2), cortn.p(3,2));
fprintf('Zeta z vs. x position:   rho = %.3f, p = %.3f\n', cortn.R(1,3), cortn.p(1,3));
fprintf('           y position:   rho = %.3f, p = %.3f\n', cortn.R(2,3), cortn.p(2,3));
fprintf('           z position:   rho = %.3f, p = %.3f\n', cortn.R(3,3), cortn.p(3,3));
fprintf('\n')
% list of successful correlations
[p(1,:),p(2,:)] = find(cortn.p<0.01);
if size(p,2)>0
    fprintf('Significant correlations to 99%% confidence:\n  ')
    coordLabels = 'xyz'; first=1;
    for i=p
        if ~first; fprintf(','); else; first=0; end
        fprintf(' zeta %s vs. %s',coordLabels(i(2)),coordLabels(i(1)));
    end
    fprintf('\n');
else
  fprintf('No significant correlations to 99%% confidence.\n');
end
fprintf('\n');

%line of best fits
lineFit = aS.lineFit;
kitLog('Results of line fitting analysis');
fprintf('\nBest fit lines for significant correlations:\n')
for i=p
  switch sign(lineFit.b(i(1),i(2)))
      case {1,0}
          sib = '+';
      case -1
          sib = '-';
  end
  fprintf('Zeta %s vs. %s:   zeta %s = %.2f%s %s %.2f, with R^2 = %.3f\n',...
      coordLabels(i(2)),coordLabels(i(1)),...
      coordLabels(i(2)),lineFit.ax(i(1),i(2)),coordLabels(i(1)),sib,abs(lineFit.b(i(1),i(2))),...
      lineFit.R2(i(1),i(2)));
end
fprintf('\n');


end %printStatistics


function aS = defaultAnalStruct

meas = struct('median',[],...
              'mean',[],...
              'stdDev',[],...
              'stdErr',[]);
stats = struct('n',0,...
               'zetaX',meas,...
               'zetaY',meas,...
               'zetaZ',meas);
correlation = struct('R',zeros(3),'p',ones(3));
lineFit = struct('ax',zeros(3),'b',zeros(3),'R2',zeros(3));

aS.version = 'ChrS analysis v1.0';
aS.date = datestr(floor(now));
aS.filename = '';
aS.filepath = '';
aS.stats = stats;
aS.correlation = correlation;
aS.lineFit = lineFit;

end %defaultAnalStruct

