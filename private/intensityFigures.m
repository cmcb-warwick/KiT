function intensityFigures(intStructs,varargin)
% Produces figures and analyses of intensity measurements for an intensity
% structure output from intensityMeasurements.
%
% Used by the kitAnalysis function.
%

% default vs user-defined options
opts.channels = [];
opts.channelLabels = {'green','red','far-red','blue'};
opts.closeFigs = 0;
opts.conditions = {'control','siRNA','DMSO','noc','taxol'};
opts.filename = [];
opts.normalise = []; % normalisation condition
opts.plotType = 'boxwhisk'; % bars, 1Dscat, boxwhisk
opts.save = 0;
opts.savePath = [];
opts = processOptions(opts,varargin{:});

% process input
if nargin<1 || isempty(intStructs)
    error('Please provide an intensity structure.')
elseif ~iscell(intStructs)
    intStructs = {intStructs};
    nConds = 1; % assume is a single condition
else
    nConds = length(intStructs);
end
% get channels analysed - NB this needs to be consistent between conditions
if isempty(opts.channels)
    chans = intStructs{1}.channels;
elseif any(~ismember(opts.channels,intStructs{1}.channels))
    error('Requested channels not all analysed in intensity structure.');
else
    chans = opts.channels;
end
% check that plotType method exists
if ~any(strcmp(opts.plotType,{'bars','1Dscat','boxwhisk'}))
    opts.plotType = 'boxwhisk';
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
    
nChans = length(chans);
chanLab = opts.channelLabels(1:nChans);
condLab = opts.conditions(1:nConds);

%% Plot figures

kitLog('Processing figures:');

% raw intensities
f = figure();
for iCond = 1:nConds
    for iChan = 1:nChans
        strCmd = ['data' num2str(iChan) '{iCond} = intStructs{iCond}.raw(:,' num2str(chans(iChan)) ');'];
        eval(strCmd);
    end
end
for iChan = 1:nChans
    strCmd = ['thisData = data' num2str(iChan) ';']; eval(strCmd);
    if ~all(isnan(cat(1,thisData{:})))
        if opts.normalise
            normFact = nanmedian(thisData{1});
            for iCond = 1:nConds
                thisData{iCond} = thisData{iCond}/normFact;
            end
        end
        subplot(1,nChans,iChan);
        switch opts.plotType
            case 'boxwhisk'
                compareBoxWhiskers(thisData(1:nConds),'withinFig',1);
            case '1Dscat'
                compare1Dscatters(thisData(1:nConds),'withinFig',1);
            case 'bars'
                compareBars(thisData(1:nConds),'withinFig',1,'errorBars','stdErr');
        end
        % aesthetics
        strTit = [chanLab{iChan} ' intensity']; title(strTit);
        if iChan==1; ylabel('background-corrected intensity');
        else; ylabel(''); end
        set(gca,'XTickLabel',condLab,'FontSize',12);
    end
end
if opts.save
  
  kitLog('Saving figure: Background-corrected intensities')
  
  % create filename
  filename = [opts.plotType '_bgCorr_' opts.filename];
  % print the figure to file
  print(f,fullfile(opts.savePath,filename),'-depsc');
  if opts.closeFigs
    close(f);
  end
  
end

% cell-normalised intensities
f = figure();
for iCond = 1:nConds
    for iChan = 1:nChans
        strCmd = ['data' num2str(iChan) '{iCond} = intStructs{iCond}.norm.cellwise(:,' num2str(chans(iChan)) ');'];
        eval(strCmd);
    end
end
for iChan = 1:nChans
    strCmd = ['thisData = data' num2str(iChan) ';']; eval(strCmd);
    if ~all(isnan(cat(1,thisData{:})))
        if opts.normalise
            normFact = nanmedian(thisData{1});
            for iCond = 1:nConds
                thisData{iCond} = thisData{iCond}/normFact;
            end
        end
        subplot(1,nChans,iChan);
        switch opts.plotType
            case 'boxwhisk'
                compareBoxWhiskers(thisData(1:nConds),'withinFig',1);
            case '1Dscat'
                compare1Dscatters(thisData(1:nConds),'withinFig',1);
            case 'bars'
                compareBars(thisData(1:nConds),'withinFig',1,'errorBars','stdErr');
        end
        % aesthetics
        strTit = [chanLab{iChan} ' intensity']; title(strTit);
        if iChan==1; ylabel('per-cell normalised intensity');
        else; ylabel(''); end
        set(gca,'XTickLabel',condLab,'FontSize',12);
    end
end
if opts.save
  
  kitLog('Saving figure: Intensities normalised to cell-by-cell control channel')
  
  % create filename
  filename = [opts.plotType '_cellNorm_' opts.filename];
  % print the figure to file
  print(f,fullfile(opts.savePath,filename),'-depsc');
  if opts.closeFigs
    close(f);
  end
  
end

% spot-normalised intensities
f = figure();
for iCond = 1:nConds
    for iChan = 1:nChans
        strCmd = ['data' num2str(iChan) '{iCond} = intStructs{iCond}.norm.spotwise(:,' num2str(chans(iChan)) ');'];
        eval(strCmd);
    end
end
for iChan = 1:nChans
    strCmd = ['thisData = data' num2str(iChan) ';']; eval(strCmd);
    if ~all(isnan(cat(1,thisData{:})))
        if opts.normalise
            normFact = nanmedian(thisData{1});
            for iCond = 1:nConds
                thisData{iCond} = thisData{iCond}/normFact;
            end
        end
        subplot(1,nChans,iChan);
        switch opts.plotType
            case 'boxwhisk'
                compareBoxWhiskers(thisData(1:nConds),'withinFig',1);
            case '1Dscat'
                compare1Dscatters(thisData(1:nConds),'withinFig',1);
            case 'bars'
                compareBars(thisData(1:nConds),'withinFig',1,'errorBars','stdErr');
        end
        % aesthetics
        strTit = [chanLab{iChan} ' intensity']; title(strTit);
        if iChan==1; ylabel('per-spot normalised intensity');
        else; ylabel(''); end
        set(gca,'XTickLabel',condLab,'FontSize',12);
    end
end
if opts.save
  
  kitLog('Saving figure: Intensities normalised to spot-by-spot control channel')
  
  % create filename
  filename = [opts.plotType '_spotNorm_' opts.filename];
  % print the figure to file
  print(f,fullfile(opts.savePath,filename),'-depsc');
  if opts.closeFigs
    close(f);
  end
  
end

%% Save results to file

if opts.save  
  % create filename for .mat struct
  filename = ['intStruct_' opts.filename];
  kitLog('Saving intensity analysis MATLAB structure: %s',filename)
  % print the intensity structure to file
  save(fullfile(opts.savePath,filename),'intStructs');
  
  kitLog('Saving .csv files:')
  for iCond = 1:nConds
    % create filename for background-corrected data
    filename = sprintf('bgcorr_%s_%s.csv',condLab{iCond},opts.filename);
    kitLog('%s',filename)
    csvwrite(fullfile(opts.savePath,filename),intStructs{iCond}.raw);
    % per-cell normalised data
    filename = sprintf('cellnorm_%s_%s.csv',condLab{iCond},opts.filename);
    kitLog('%s',filename)
    csvwrite(fullfile(opts.savePath,filename),intStructs{iCond}.norm.cellwise);
    % per-spot normalised data
    filename = sprintf('spotnorm_%s_%s.csv',condLab{iCond},opts.filename);
    kitLog('%s',filename)
    csvwrite(fullfile(opts.savePath,filename),intStructs{iCond}.norm.spotwise);
  end
end


end