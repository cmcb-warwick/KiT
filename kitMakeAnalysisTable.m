function analysisTable = kitMakeAnalysisTable(jobset,saveToCSV, channel, dt)
%% kitMakeAnalysisTable(jobset,saveToCSV, channel)
%%
%% jobset - processed paired jobset object
%% saveToCSV - logical, should output be saved
%% channel - int, relevant channel
%%
%%Take a processed, paired jobset file and convert analysis to a convenient
%%table format
%% Will account for filtered spots and ignore any that have been filtered out
%%Jonathan U Harrison 2019-04-12
%%%%%%%%%%%%%

if nargin < 2
    saveToCSV = 0; %save output as csv format
end

if nargin < 3
    channel = 1;
end

if nargin < 4
    dt = [];
end

%get the data
job = kitLoadAllJobs(jobset);
nMovs = length(job);

for jobInd = 1:nMovs
    try
    dataStruct = job{jobInd}.dataStruct{channel};
    
    if ~isfield(job{jobInd},'dataStruct')
        error('Expect a jobset file that has already been processed and paired');
    end
    
    if ~isfield(dataStruct,'failed') || dataStruct.failed
        warning('Movie %d has failed so skipping \n',jobInd);
        continue
    end
    
    nFrames = job{jobInd}.metadata.nFrames;
    if nFrames > 1 %treat time sereies and single time points differently
        varnames = {'Position','Amplitude','Frame','Time','SisterPairID','SisterID'};
        analysisTable = cell2table(cell(0,7),'VariableNames',cat(2,varnames,'movieID'));
    else
        varnames = {'Position','Amplitude'};
        if isfield(job{1},'categories')
            numVars = 3 + numel(fieldnames(job{1}.categories));
            analysisTable = cell2table(cell(0,numVars),'VariableNames',cat(2,varnames,fieldnames(job{1}.categories)),'movieID');
        else
            numVars = 3;
            analysisTable = cell2table(cell(0,numVars),'VariableNames',cat(2,varnames,'movieID'));
        end
    end
    if nFrames > 1
        warning('Categories cannot yet be used for movies. Making a basic table')
        if isfield(dataStruct,'sisterList')
            nSisters = length(dataStruct.sisterList);
            %create columns for table
            position = zeros(nFrames*nSisters*2,3);
            amplitude = zeros(nFrames*nSisters*2,3);
            frame = repmat(repmat((1:nFrames)',nSisters,1),2,1);
            sisterPairID = repmat(1:nSisters,nFrames,1);
            sisterPairID = repmat(sisterPairID(:),2,1);
            sisterID = [ones(nFrames*nSisters,1);2*ones(nFrames*nSisters,1)];
            for k = [1,2]
                for j = 1:nSisters
                    for i = 1:nFrames
                        if k==1 %use coords1 or coords2
                            position((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                                dataStruct.sisterList(j).coords1(i,1:3);
                            amplitude((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                                dataStruct.sisterList(j).coords1(i,4:6);
                        elseif k==2
                            position((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                                dataStruct.sisterList(j).coords2(i,1:3);
                            amplitude((k-1)*nFrames*nSisters + (j-1)*nFrames + i,:) = ...
                                dataStruct.sisterList(j).coords2(i,4:6);
                        end
                    end
                end
            end
            if isempty(dt)
              dt = dataStruct.dataProperties.frameTime(1,2)-dataStruct.dataProperties.frameTime(1,1);
            end
            T = table(position,amplitude,frame,(frame-1)*dt,sisterPairID,sisterID,...
                'VariableNames', varnames);
        else
            error('No sisterList found but need paired joblist to make table from');
        end
    else
        %put all spots into a table structure
        T = table(dataStruct.initCoord.allCoord(:,1:3), ...
            dataStruct.initCoord.allCoord(:,4:6), ...
            'VariableNames',{'Position','Amplitude'});
        % add category information from manual labelling via gui
        if isfield(job,'categories') %have added some kind of categories
            categoryNames = fieldnames(job.categories);
            for j=1:numel(fieldnames(job.categories))
                T.(categoryNames{1}) = zeros(size(T,1),1);
                T.(categoryNames{1})(job.categories.(categoryNames{1})) = 1;
            end
        end
    end
    T.movieID = jobInd*ones(size(T,1),1);
    analysisTable = [analysisTable; T];
    if saveToCSV
        outname = kitGenerateOutputFilename(job{jobInd});
        outname = strcat(outname(1:(end-3)),'csv'); %replace file ending of mat with csv
        writetable(analysisTable, outname);
    end
    catch
        fprintf('Unsuccessful in converting to table: job %d run correctly?\n',jobInd);
        continue
    end
    fprintf('Successful in converting to table for job %d.\n',jobInd);
end

