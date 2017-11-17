function printHeader(head,fid)

if nargin<2
  fid=0;
end

under = repmat('-',[1 length(head)]);

fprintf(fid,'%s\n',head);
fprintf(fid,'%s\n',under);
