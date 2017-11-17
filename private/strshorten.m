function s=strshorten(s,l)
% Shorten string to length l using dots in middle.

if length(s) > l
  m = floor(l/2) - 2;
  n = m + mod(l,2);
  s = [s(1:m) '...' s(end-n:end)];
end
