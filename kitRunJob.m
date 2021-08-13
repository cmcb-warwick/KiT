function kitRunJob(jobset,varargin)
% KITRUNJOB Runs tracking analysis on jobset
%
%    KITRUNJOB(JOBSET,...) Runs tracking analysis on JOBSET. JOBSET should be
%    output from either KITSETUPJOBS or KITLOADJOBSET. Additional options can
%    be supplied as string/value pairs.
%
%    Options (default in {}):-
%
%    tasks: {1:9} or a vector. Specify tasks to perform on jobs:
%                              1: finding spots
%                              2: fitting plane
%                              3: tracking spots
%                              4: grouping sisters
%                              5: extracting tracks
%                              6: neighbour spot finding
%                              7: updating classes
%                              8: aligning
%                              9: intensity
%
%    existing: {0} or 1. Load existing jobs first.
%
%    exec: {'serial'}, 'batch', 'pbs', 'pbsstage'. Execution mode.
%
%    subset: {[]} or vector of job numbers. Set to a subset of indices of movies
%    to analysis, instead of processing them all.
%
%    errorfail: {0} or 1. Crash out if error occurs if 1.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

% Check minimum MATLAB version.
% FIXME Check minimum toolbox versions also.
if verLessThan('matlab','8.6')
  error('Minimum required MATLAB version is 8.6 (R2015b)');
end

% Download BioFormats, if required.
kitDownloadBioFormats();

% Upgrade jobset options, if required.
if ~isfield(jobset,'jobsetVersion') || ...
    jobset.jobsetVersion < kitVersion(2)
  jobset = kitJobset(jobset);
end

nROIs = length(jobset.ROI);

% Default options.
options.subset = 1:nROIs;
options.exec = 'serial';
options.errorfail = 0;
options.tasks = 1:9;
options.existing = 0;
options.callback = [];
options.email = [];
% Get user options.
options = processOptions(options, varargin{:});

% Check options.
if ~all(ismember(options.subset,1:nROIs))
  error('Subset values must be in range 1 to %d',nROIs);
end
switch jobset.options.jobProcess
    case 'zonly'
        options.tasks = setdiff(options.tasks,[3,4,5,7,8]);
    case 'chrshift'
        options.tasks = setdiff(options.tasks,[2,3,4,5,7,8]);
end
if all(~jobset.options.intensity.execute)
  options.tasks = setdiff(options.tasks,9);
end
if ~ismember(options.tasks,1)
    options.existing = 1;
end

% If using matlabpool for parallel computation, report workers.
[~,name] = fileparts(jobset.filename);
switch options.exec
  case 'batch'
    kitLog(['Running ' name ' in parallel']);
  case 'serial'
    kitLog(['Running ' name ' serially']);
  case {'pbs','pbsstage'}
    kitLog(['Running ' name ' using PBS']);
end

if isfield(jobset,'variantName')
  fprintf('Jobset variant: %s\n',jobset.variantName);
end

% Copy out job info for each movie.
jobs = cell(nROIs,1);
for i= options.subset
  if options.existing
    jobs{i} = kitLoadJob(jobset,i);
    % Copy over any new options.
    jobs{i}.options = jobset.options;
    if length(jobs{i}.ROI)>1
        jobs{i}.ROI = jobset.ROI(i);
    end
  else
    jobs{i} = jobset;
    jobs{i}.index = i;
    jobs{i}.ROI = jobset.ROI(i);
  end
  % Update versions, may be different to jobset creator.
  jobs{i}.version = kitVersion();
  jobs{i}.matlabVersion = version;
  %jobs{i}.lociVersion = char(loci.formats.FormatTools.VERSION);
  % Record host.
  if ispc
    [~,jobs{i}.host] = system('echo %COMPUTERNAME%');
  else
    [~,jobs{i}.host] = system('hostname');
  end
end


exceptions = [];
if strcmp(options.exec,'pbs')
  cmd = sprintf('qsub -N KiT_%s -v JOBSET_FILE="%s" -t %s',name,jobset.filename,strjoin(arrayfun(@num2str,options.subset,'uniformoutput',0),','));
  if ~isempty(options.email)
    cmd = [cmd ' -m e -M ' options.email];
  end
  cmd = [cmd ' private/pbstemplate.pbs'],
  [status,result] = system(cmd);
  if status~=0
    error('Error submitting PBS job: %s',result);
  end
else
  for i = options.subset
    switch options.exec
      case 'batch'
        kitLog('Submitting tracking job %d', i);
        batchJob{i} = batch(@kitTrackMovie, 1, {jobs{i},options.tasks});
      case 'serial'
        try
          kitLog('Tracking job %d', i);
          kitTrackMovie(jobs{i},options.tasks);
          if ~isempty(options.callback)
            options.callback(i);
          end
        catch me
          kitLog('Error in job %d: %s',i,me.identifier);
          ex.me = me;
          ex.idx = i;
          exceptions = [exceptions ex];
          if options.errorfail
            disp(getReport(me));
            throw(me);
          end
        end
      case 'pbsstage'
        kitLog('Submitting tracking job %d to PBS with staging',i);
        [~,trackFile,ext] = fileparts(kitGenerateOutputFilename(jobs{i}));
        trackFile = [trackFile ext];
        cmd = sprintf('qsub -N KiT_%s_%d -v MOVIE_FILE="%s",MOVIE_DIR="%s",TRACK_FILE="%s",JOBSET_FILE="%s",JOB_ID=%d private/pbstemplatestage.pbs',name,i,jobset.ROI(i).movie,jobset.movieDirectory,trackFile,jobset.filename,i);
        [status,result] = system(cmd);
        if status~=0
          error('Error submitting PBS job for job %d: %s',i,result);
        end
        pause(0.5); % Wait a little bit.
    end
  end
end


if ~isempty(exceptions)
  disp('Errors occured:')
end
for i = 1:length(exceptions)
  ex = exceptions(i);
  fprintf('In job %d, error %s:\n',ex.idx,ex.me.identifier);
  disp(getReport(ex.me));
end

if strcmp(options.exec,'batch')
  % Wait for parallel tracking jobs.
  for i = options.subset
    kitLog('Waiting for %d tracking jobs to complete', length(options.subset)-i+1);
    b = batchJob{i};
    wait(b);
    diary(b);
    delete(b);
  end
end

switch options.exec
  case {'serial','batch'}
    kitLog('Tracking complete');
  case 'pbs'
    kitLog('Submission complete');
end

% Dump jobset diagnostics.
[pathstr,name,ext] = fileparts(jobset.filename);
diagfile = ['diags_' name];
if isfield(jobset,'variantName')
  diagfile = [diagfile '_' jobset.variantName];
end
diagfile = fullfile(pathstr,[diagfile '.txt']);
fid = fopen(diagfile,'wt');
C = onCleanup(@() fclose(fid));
c = 1:length(jobset.options.coordMode);
cindx = ~strcmp(jobset.options.coordMode, 'none');
kitJobsetDiagnostics(jobset,c(cindx),0,fid);
