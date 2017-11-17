function kitDownloadBioFormats()
% KITDOWNLOADBIOFORMATS Download BioFormats, if required.
%
%     Download BioFormats, if required. Requires current directory to be KiT
%     directory.
%
% Copyright (c) 2015 Jonathan W. Armond

if exist('bfGetPlane') == 2
  return % already available
end

filename = 'bfmatlab';
url = ['http://downloads.openmicroscopy.org/bio-formats/5.1.7/artifacts/' ...
       'bfmatlab.zip'];

% Check for existing distribution.
download = 0;
if exist(filename,'dir') == 0
  % Directory doesn't exist, check for zip.

  if exist([filename '.zip'],'file') == 0
    % Zip doesn't exist. Download it.
    kitLog('Downloading BioFormats...');
    % Check we are in KiT directory.
    if exist(fullfile(pwd,'kitDownloadBioFormats.m'),'file') ~= 2
      error('Must be in KiT directory to download BioFormats');
    end

    unzip(url);
    % Check successful.
    addpath bfmatlab;
    if bfCheckJavaPath(1)
      kitLog('Download complete.');
    else
      error('Failed to download BioFormats. Please download bfmatlab.zip from http://downloads.openmicroscopy.org/bio-formats/ and unzip in the KiT directory.');
    end
  else
    % Unzip.
    unzip([filename '.zip']);
  end
end
