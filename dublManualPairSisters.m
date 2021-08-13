function dublManualPairSisters(jobset,varargin)
% DUBLMANUALPAIRSISTERS Manually select pairs of kinetochores.
%
%    JOB = DUBLMANUALPAIRSISTERS(JOB,...) Images of each kinetochore found
%    in a JOB are shown and their centre plotted (in white). Centres of
%    nearby kinetochores (potential sisters, in green) are also shown, so
%    that the user can click next to the spot they believe to be the
%    original spot's sister.
%    Kinetochores with no clear sister can be omitted from the sister list
%    by pressing any key.
%    Nearby kinetochores already allocated (red), and kinetochores
%    previously decided to have no clear pair (yellow), are also plotted.
%
%    Options, defaults in {}:-
%
%    coordChans: {[1 2]} or any subset of [1,2,3]. The channels for which
%                   sisterList structures will be produced.
%
%    imageChans: {[1 2]} or any subset of [1,2,3]. The channels which will
%                   be plotted when inspecting spots for pairing.
%
%    maxSisSep: {2.5} or distance in µm. Maximum distance between the
%                   original spot and potential sisters.
%
%    metaphasePlate: {0} or 1. Whether or not to plot a metaphase plate.
%                   N.B. This doesn't work yet. Future plan.
%
%    mode: {'dual'}, 'zoom' or 'full'. Whether to show images as:
%                   - 'full'
%                     full images.
%                   - 'zoom'
%                     zooming in on kinetochores to an area defined by the
%                     'maxSisSep' option.
%                   - 'dual'
%                     a combination of the two.
%
%    plotChan: {1}, 2 or 3. The channel which will be used for plotting
%                   spot centres when inspecting spots for pairing. This
%                   will be stored in .options.
%                   N.B. This only works for one channel so far, but may
%                   be increased in future.
%
%    redo: {0} or 1. Whether to completely redo sister pairing, i.e. don't
%                   plot spots as being previously allocated etc.
%
%    subset: {[]} or list of movie numbers. Sub-list of movies to process.
%
%    verbose: {0} or 1. Whether or not to print progress to the command
%                   line.
%
%    zProjRange: {2} or distance in pixels. Total distance about the
%                   original spot over which to project in the
%                   z-coordinate.
%
%
%    Future directions:-
%
%    - Incorporate plotting of metaphase plate
%    - Allow processing of movies as well as single z-stacks
%       - Pair tracks instead of spots, show dragon tails of nearby tracks
%       - Start at centre of track being paired
%       - Allow two tracks to be paired to another track, but no overlaps
%         allowed
%
% Copyright (c) 2016 C. A. Smith

% default options
opts.chanOrder = [2 1 3]; % order of the channels in the movie
opts.contrast = {[0.1 0.9995] [0.1 0.9995] [0.1 1]};
opts.imageChans = 1; % imaging more channels may provide more information
opts.maxSisSep = 2.5; % maximum distance in µm between potential sisters being plotted
opts.mode = 'dual';
opts.plotChan = 1; % channel coordinates to be plotted for allocation
opts.redo = 0;
opts.subset = [];
opts.verbose = 0;
opts.zProjRange = 2; % number of pixels over which to project in z

% process user-defined options
opts = processOptions(opts,varargin{:});

% turn off all warnings
warning('off','all')

%% FRONT END

% get jobs
jobs = kitLoadAllJobs(jobset);

% if no subset given, use full length of jobs
if isempty(opts.subset)
    opts.subset = 1:length(jobs);
end

% start counter
counter = 1;
nMovies = length(opts.subset);
userStatus = 'completed';

% loop over movies to run kitManualPairSisters
for iMov = opts.subset
    
    % output log information
    kitLog('Pairing kinetochores in movie %i (%i of %i):',iMov,counter,nMovies);
    counter = counter+1;
    
    % find whether is a single z-stack, or a time series
    % perform sister pairing
%     if jobs{iMov}.metadata.nFrames > 1
%       [jobs{iMov},userStatus] = pairTracks(jobs{iMov},opts);
%     else
      [jobs{iMov},userStatus] = pairSpots(jobs{iMov},opts);
%     end
    
    if strcmp(userStatus,'userPaused')
        [~,lastMovie] = find(opts.subset == iMov);
        if lastMovie == 1
            lastMovie = 0;
            break
        else
            lastMovie = opts.subset(lastMovie-1);
            break
        end
    end
end

% check user input
switch userStatus
    case 'userPaused'
        if lastMovie == 0
            kitLog('User paused sister pairing before any movies saved.')
            return
        else
            kitLog('User paused sister pairing. Sister pairing saved up to movie %i.',lastMovie)
        end
    case 'completed'
        kitLog('Sister pairing completed for all movies.')
end

end
