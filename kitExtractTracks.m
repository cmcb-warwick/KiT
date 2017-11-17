function job=kitExtractTracks(job, channel)
%KITEXTRACTTRACKS Extracts individual tracks from dataStructs.
%
% SYNOPSIS: job=kitExtractTracks(job, channel)
%
% INPUT job: Struct containing tracking job results.
%       channel: Channel number to operate on.
%
% OUTPUT job: As input but with new field:
%                 .trackList
%
% Copyright (c) 2012 Jonathan W. Armond

dataStruct = job.dataStruct{channel};
tracks = dataStruct.tracks;
planeFit = dataStruct.planeFit;
sisterList = dataStruct.sisterList;
nTracks = length(tracks);
trackList(1:nTracks) = struct(...
    'coords',nan(job.metadata.nFrames,6),...
    'featIndx',nan(job.metadata.nFrames,1),...
    'attach',zeros(job.metadata.nFrames,1));
opts = job.options;

% Extract track coordinates.
for j=1:nTracks
    track = tracks(j);
    featIndx = track.tracksFeatIndxCG;
    startTime = track.seqOfEvents(1,1);
    endTime = track.seqOfEvents(2,1);

    for t=startTime:endTime
        featT = t-startTime+1;
        if featIndx(featT) ~= 0
            trackList(j).coords(t,:) = planeFit(t).rotatedCoord(featIndx(featT),:);
            trackList(j).featIndx(t) = featIndx(featT);
        end
    end

    % Assign P or AP directions.
    trackList(j).direction = kitAssignDirection(trackList(j).coords(:,1),...
                                                'exp',job.options.dirAssignExpWeight);

    % Set sister to empty in case no sisters.
    trackList(j).sister = [];
end

% Use sisters to assign poles. Standardize so that +ve X displacement is
% towards pole.
trackPairs = sisterList(1).trackPairs;
if isempty(sisterList(1).trackPairs)
  nSisters = 0;
else
  nSisters = length(sisterList);
end
done = [];
for i=1:nSisters
    a = trackPairs(i,1);
    b = trackPairs(i,2);
    if a > nTracks || b > nTracks
        warning('Sister tracks out of range');
        continue;
    end

    % Find mean x pos of each track.
    aPos = nanmean(trackList(a).coords(:,1));
    bPos = nanmean(trackList(b).coords(:,1));

    % + => attached to x>0 pole, - => attached to x<0 pole.
    % 2 => pole assigned via sister pairing.
    % 1 => pole assigned based on track extent.
    if aPos < bPos
        trackList(a).attach = -2;
        trackList(a).direction = -trackList(a).direction;
        trackList(b).attach = +2;
    else
        trackList(b).attach = -2;
        trackList(b).direction = -trackList(b).direction;
        trackList(a).attach = +2;
    end
    trackList(a).sister = b;
    trackList(b).sister = a;

    % Mark tracks as resolved.
    done = [done; a; b];
end

if opts.debug.asserts
  % Check trackList with sisters matches sisterList
  for i=1:size(trackPairs,1)
    a = trackPairs(i,1);
    b = trackPairs(i,2);
    notnan = ~isnan(sisterList(i).coords1(:,1));
    if ~isequal(trackList(a).coords(notnan,:),sisterList(i).coords1(notnan,:)) || ...
        ~isequal(trackList(b).coords(notnan,:),sisterList(i).coords2(notnan,:))
      error('trackList(%d,%d) <-> sisterList(%d) mismatch',a,b,i);
    end
  end
end

% For remaining tracks assess attachment based on max or min pos.
remain = setdiff(1:nTracks,done);
for i=1:length(remain)
    xCoord = trackList(remain(i)).coords(:,1);
    [~,minIdx] = nanmin(abs(xCoord));
    [~,maxIdx] = nanmax(abs(xCoord));
    minPos = xCoord(minIdx);
    maxPos = xCoord(maxIdx);

    if minPos > maxPos
      trackList(remain(i)).attach = -1;
      trackList(remain(i)).direction = -trackList(remain(i)).direction;
    else
      trackList(remain(i)).attach = +1;
    end
end

if opts.debug.asserts
  % Check if sisters attached to -ve pole are nearer to it than other sister.
  for i=1:nTracks
    t = trackList(i);
    if isempty(t.sister)
      continue;
    end
    s = trackList(t.sister);
    if sign(t.attach) ~= sign(nanmean(t.coords(:,1))-nanmean(s.coords(:,1)))
      error('sister pole attachment inconsistent with trajectory (%d,%d)',i,t.sister);
    end
  end

  % Check P displacements are x +ve and AP are x -ve for x>0 attached KT, and
  % vice versa.
  % for i=1:nTracks
  %   t = trackList(i);
  %   if sign(t.attach)*median(diff(t.coords(t.direction==1,1))) < 0
  %     error('P displacements are opposite sign (%d)',i);
  %   end
  %   if sign(t.attach)*median(diff(t.coords(t.direction==-1,1))) > 0
  %     error('AP displacements are opposite sign (%d)',i);
  %   end
  % end
end

% Record changes.
job.dataStruct{channel}.trackList = trackList;
