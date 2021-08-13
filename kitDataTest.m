% Download test data to evaluate algorithms on
% Note this will not be able to be automated unless movie included ...
% with repo or code altered
%
% Jonathan U Harrison 2019-02-06
%%%%%%%%%%
%First check if Metaphase example is present in test directory, and if not
%then download
checkDirCode = exist('testdata','dir');
if checkDirCode==7
    fprintf('\ntest directory exists, checking for Metaphase example movie\n');
elseif checkDirCode==0
    fprintf('\nDirectory does not exist. Attempting to create. \n');
    success = mkdir('testdata');
    if ~success
        error('Directory creation seemed to go wrong');
    end
else
    error('Cannot find /testdata. Instead conflicting test object found. Check working directory and/or workspace?');
end    

if exist('testdata/MetaphaseExample.ome.tiff','file')
    fprintf('\nMetaphase example exists. using this to perform further tests\n');
else
    fprintf('\nFirst downloading Metaphase movie to use for further tests\n'); 
    fprintf('\nPlease save the file in the testdata folder\n');     
    web('https://bitbucket.org/jarmond/kit/downloads/MetaphaseExample.ome.tiff','-browser')
    pause(120);
    success = movefile('~/Downloads/MetaphaseExample.ome.tiff',...
        'testdata/MetaphaseExample.ome.tiff');   
    fprintf('\nDownload successful \n');
    if ~exist('testdata/MetaphaseExample.ome.tiff','file') || ~success
        error('\n Downloading did not seem to work this time. Check internet connection?'); 
    end
end
%%%%%%%%%
%bioformats?

% try opening the file
[metadata, reader] = kitOpenMovie('testdata/MetaphaseExample.ome.tiff');

% Then come the tests
%%%%%%%%%%%%%%
%% Test1: is reader object empty?
assert(~isempty(reader),'Hopefully we read something');

%% Test2: check metadata
assert(isstruct(metadata),'metadata exists and is a struct');
assert(isfield(metadata,'nFrames'));
assert(isfield(metadata,'nChannels'));
assert(metadata.nFrames==150,'check number of frames');
assert(metadata.nChannels==1,'check number of channels');




