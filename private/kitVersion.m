function version=kitVersion(mode)
% KITVERSION Return string containing version number

if nargin < 1
  mode = 1;
end

switch mode
  case 1
    % KiT version.
    version = '2.1.11';
  case 2
    % Jobset structure version.
    version = 8;
  otherwise
    error('Unknown version mode');
end
