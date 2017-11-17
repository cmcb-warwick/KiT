function allJobs=kitNormalizeJobs(jobs,normalize)
%KITCOLLATEJOBS Combines multiple jobs into one and normalizes intensities.
%
% SYNOPSIS: allJobs=kitNormalizeJobs(jobs)
%
% INPUT jobs: Cell array of jobs.
%
%       normalize: 1 do background subtraction and scale, 2 just background subtraction
%
% OUTPUT allJobs: Struct similar to job structs but containing combined data for
%                 for the following fields:
%                 .sisterList
%                 .trackList
%                 .trackInt
%
% Copyright (c) 2014 Jonathan W. Armond

if nargin<2
  normalize=1;
end

n = length(jobs);
cellCtr = 1;
trackCtr = ones(3,1); % assume max 3 channels
sisterCtr = trackCtr;
fields = {'intensity','intensity_median','intensity_min','intensity_max','intensity_ratio'};

for i=1:n

  job = jobs{i};
  numCh = length(job.dataStruct);
  for c = 1:numCh
    if ~isempty(job.dataStruct{c}) && ...
          isfield(job.dataStruct{c},'cellInt') && ~isempty(job.dataStruct{c}.cellInt.back)
        % Background intensity is stored in cellInt.
        %back = mean(job.dataStruct{c}.cellInt.backMode,1);
        %backDiff = mean(job.dataStruct{c}.cellInt.backDiff,1);
        back = job.dataStruct{c}.cellInt.backMode(1,:);
        backDiff = job.dataStruct{c}.cellInt.backDiff(1,:);
        allJobs.dataStruct{c}.background(cellCtr,:) = back;
        allJobs.dataStruct{c}.backgroundDiff(cellCtr,:) = backDiff;

        % Per cell data
        allJobs.dataStruct{c}.cellInt(cellCtr) = job.dataStruct{c}.cellInt;
        allJobs.dataStruct{c}.movie{cellCtr} = job.movie;

        % Augment track pair indices.
        if ~isempty(job.dataStruct{c}.sisterList(1).trackPairs)
          job.dataStruct{c}.sisterList(1).trackPairs(:,1:2) = ...
              job.dataStruct{c}.sisterList(1).trackPairs(:,1:2) + ...
              (trackCtr(c)-1);
        end

        for j = 1:length(job.dataStruct{c}.trackList)
          k = trackCtr(c);
          % Augment track sister index.
          %job.dataStruct{c}.trackList(j).sister = ...
          %    job.dataStruct{c}.trackList(j).sister + (trackCtr(c)-1);
          % Normalize intensity.
          trackInt = job.dataStruct{c}.trackInt(j);
          if normalize
            for f = 1:length(fields)
              for ii = 1:size(trackInt.intensity,2)
                if normalize == 2
                  trackInt.(fields{f})(:,ii) = (trackInt.(fields{f})(:,ii) - ...
                                              back(ii));
                else
                  trackInt.(fields{f})(:,ii) = (trackInt.(fields{f})(:,ii) - ...
                                              back(ii)) / backDiff(ii);
                end
              end
            end
          end

          allJobs.dataStruct{c}.trackInt(k) = trackInt;
          allJobs.dataStruct{c}.trackList(k) = job.dataStruct{c}.trackList(j);
          allJobs.dataStruct{c}.jobIdx(k) = i;
          allJobs.dataStruct{c}.trackIdx(k) = j;
          trackCtr(c) = trackCtr(c)+1;
        end
        for j = 1:length(job.dataStruct{c}.sisterList)
          k = sisterCtr(c);
          allJobs.dataStruct{c}.sisterList(k) = job.dataStruct{c}.sisterList(j);
          allJobs.dataStruct{c}.sisterIdx(k) = j;
          sisterCtr(c) = sisterCtr(c)+1;
        end
        cellCtr = cellCtr+1;
    end
  end
end
