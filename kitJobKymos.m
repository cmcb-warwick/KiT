function kitJobKymos(job,channel,outdir,varargin)
% KITJOBKYMOS generates kymographs for all sister pairs in job

mkdir(outdir);

ds = job.dataStruct{channel};
if isempty(ds.sisterList(1).trackPairs)
  disp('No sisters.');
  return
end
for i=1:length(ds.sisterList)
  tracks = ds.sisterList(1).trackPairs(i,1:2);
  fprintf('Generating kymograph for sister %d (tracks %d,%d)\n',i,tracks(1), ...
          tracks(2));
  outfile = fullfile(outdir,sprintf('kymo-sp%02d.png',i));
  kitMakeKymograph(job,i,'trackChannel',channel,'outfile',outfile,varargin{:});
end
