function kitLog(fmt, varargin)
% KITLOG Print formatted logging message
%
%    KITLOG(FMT,...) Print formatted logging message. FMT is a printf-style
%    format string. Supply appropriate number of arguments after FMT.
%
% Copyright (c) 2012 Jonathan W. Armond

% Add timestamp to message.
fprintf([datestr(now, 'HH:MM:SS') ': ' fmt '\n'], varargin{:});
