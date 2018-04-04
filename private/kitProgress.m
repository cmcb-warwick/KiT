function lastUpdate = kitProgress(progress, lastUpdate, updateInterval)
% KITPROGRESS Print formatted progress message
%
% SYNOPSIS: kitProgress(progress, lastUpdate, updateInterval)
%
% INPUT progress: Fraction complete [0,1].
%
%       lastUpdate: Last timestamp progress was called with (optional).
%                   Omit for first call.
%
% Modified by: Christopher A. Smith
% Copyright (c) 2012 Jonathan W. Armond

currentClock = clock;

if nargin<3
    updateInterval = 10.0; % seconds
end

% Add timestamp to message.
fmt = [datestr(now, 'HH:MM:SS') ': %6.2f%% complete\n'];

msg = sprintf(fmt, 100*progress);
msg = deblank(msg);
if progress == 1
    % Print message.
    fprintf(repmat('\b',1,length(msg)+1));
    fprintf('%s\n',msg);
elseif nargin < 2 || etime(currentClock, lastUpdate) > updateInterval
    % Print message.
    if progress > 0
      fprintf(repmat('\b',1,length(msg)+1));
    end
    fprintf('%s\n',msg);
end
lastUpdate = currentClock;