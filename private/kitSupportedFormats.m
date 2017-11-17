function fmtregexp=kitSupportedFormats(asfilterspec)
% KITSUPPORTEDFORMATS Regular expression matching all supported formats.
%
%    FMTREGEXP = KITSUPPORTEDFORMATS() Returns regular expression matching all
%    supported formats.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<1
  asfilterspec = 0;
end

% Use bfGetFileExtensions
addpath bfmatlab;
bfCheckJavaPath(1);
exts = bfGetFileExtensions;
if asfilterspec == 0
  allExt = exts{1,1};
  allExt = strrep(allExt,';','|');
  allExt = strrep(allExt,'*','');
  allExt = strrep(allExt,'.','\.');
  fmtregexp = ['(' allExt ')$'];
else
  fmtregexp = exts;
end
