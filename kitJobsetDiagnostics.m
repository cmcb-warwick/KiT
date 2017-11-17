function kitJobsetDiagnostics(jobset,channel,short,fid)
% KITJOBSETDIAGNOSTICS Prints a summary of tracking results for all jobs in jobset.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<2
  channel = 1;
end
if nargin<3
  short = 1;
end
if nargin<4
  fid = 1;
end

id=@(x) x; % identity

f = struct('field','job','fn',id,'name','job','fmt','d');
f(end+1) = struct('field','nSpotsPerFrame','fn',id,'name','# spots','fmt','.1f');
f(end+1) = struct('field','nTracks','fn',id,'name','# tracks','fmt','d');
f(end+1) = struct('field','nSisters','fn',id,'name','# sisters','fmt','d');
f(end+1) = struct('field','avgSisterTrackLength','fn',id,'name','sister length','fmt','.1f');
f(end+1) = struct('field','sisterPoints','fn',id,'name','sister points','fmt','d');
f(end+1) = struct('field','percentWithPlane','fn',id,'name','plane fits','fmt','.1f');
f(end+1) = struct('field','sisterVar','fn',@(x) sqrt(x),'name','std dev','fmt','.3f');
f(end+1) = struct('field','sisterDisp','fn',id,'name','mean disp','fmt','.3f');
f(end+1) = struct('field','nLongSisters','fn',id,'name','# long','fmt','d');
f(end+1) = struct('field','nFullSisters','fn',id,'name','# full','fmt','d');
f(end+1) = struct('field','elapsedTime','fn',@(x) floor(x/60),'name','time(min)','fmt','d');
fields = f;

% Read in diagnostics from each job.
nFields = length(fields);
nJobs = length(jobset.ROI);
stats = nan(nJobs,nFields);
for i=1:nJobs
  try
    job = kitLoadJob(jobset,i);
  catch
    continue;
  end
  if ~isfield(job,'dataStruct')
    continue
  end

  ds = job.dataStruct{channel};
  if isempty(ds)
    continue
  end

  if ~isfield(ds,'diagnostics')
    continue
  end
  diag = ds.diagnostics;
  stats(i,1) = i;
  for f=2:nFields
    field = fields(f).field;
    if isfield(diag,field)
      stats(i,f) = fields(f).fn(diag.(field));
    else
      warning('Missing field: %s',field);
    end
  end
end

% Print out diagnostics table.
printHeader(strjoin({fields.name},'|'),fid);
l = cellfun(@(x) length(x), {fields.name});
for i=1:size(stats,1)
  for j=1:size(stats,2)
    fmt = sprintf('%%%d%s ',l(j),fields(j).fmt);
    fprintf(fid,fmt,stats(i,j));
  end
  fprintf(fid,'\n');
end

% Print totals
fprintf(fid,'sum ');
for j=2:size(stats,2)
  fmt = sprintf('%%%d.1f ',l(j));
  fprintf(fid,fmt,nansum(stats(:,j)));
end
fprintf(fid,'\n');

% Print averages
fprintf(fid,'avg ');
for j=2:size(stats,2)
  fmt = sprintf('%%%d.1f ',l(j));
  fprintf(fid,fmt,nanmean(stats(:,j)));
end
fprintf(fid,'\n');
