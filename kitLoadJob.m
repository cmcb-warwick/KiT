function job=kitLoadJob(jobset, index)
% KITLOADJOB Load tracking data from job
%
%    JOB = KITLOADJOB(JOBSET,INDEX) Load tracking data from job number INDEX
%    from JOBSET.
%
% Copyright (c) Jonathan Armond 2013

% Generate output name.
job = jobset;
job.movie = jobset.ROI(index).movie;
job.index = index;
outputName = kitGenerateOutputFilename(job);

% Load file.
job = load(outputName);

% Update paths to match jobset. Data may have moved.
job.movieDirectory = jobset.movieDirectory;
job.filename = jobset.filename;
job.output = outputName;
