function kitPrintDiagnostics(dataStruct,process)
% KITPRINTDIAGNOSTICS Print out basic tracking diagnostics
%
% Created by: Jonathan W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

diag = dataStruct.diagnostics;

fprintf('Elapsed compute time: %s\n',timeString(diag.elapsedTime));
if diag.failed
  fprintf('Job failed');
else
  fprintf('Particles per frame: %.1f\n',diag.nSpotsPerFrame);
  if strcmp(process,'zandt')
    fprintf('# sister pairs tracked: %d\n',diag.nSisters);
    fprintf('# individual KTs tracked: %d\n',diag.nTracks);
    fprintf('Average sister pair track length: %.1f\n',diag.avgSisterTrackLength);
    fprintf('Average KT track length: %.1f\n',diag.avgTrackLength);
    fprintf('# long sisters (75%% length): %d\n',diag.nLongSisters);
    fprintf('# full sisters (100%% length): %d\n',diag.nFullSisters);
  end
  if ~strcmp(process,'chrshift')
    fprintf('Frames with a plane fit: %.2f%%\n',diag.percentWithPlane);
  end
end