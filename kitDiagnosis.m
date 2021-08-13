%wraps jonthanDiagnoseAndTest
%should go through given folder to find kittracking files and run diagnosis
%for these
%Jonathan U Harrison 2019-02-15
%%%%%%%%%%%%

fileList = dir('../../Data/2019-04-23/kittracking*v101*capture9*.mat');
% folder = uigetdir();
% fileList = dir(fullfile(folder, '*kittracking*.mat'));
nFiles = size(fileList,1);
fprintf('Found %d files. \n', nFiles);
for j=1:nFiles
    pathToFile = sprintf('%s/%s',fileList(j).folder,fileList(j).name);
    fprintf('Processing %dth file: %s\n',j,pathToFile);
    job = kitLoadJob(pathToFile,1);
    jonathanDiagnoseAndTest(job,1,0);
    fprintf('Done \n');
end

