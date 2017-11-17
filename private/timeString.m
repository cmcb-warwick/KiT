function str=timeString(s)
% Convert seconds to string of days, hours, mins, seconds.

d = 86400;
h = 3600;
m = 60;

days = floor(s/d);
s = rem(s,d);
hours = floor(s/h);
s = rem(s,h);
mins = floor(s/m);
secs = rem(s,m);

str = [];
if days > 0
  str = sprintf('%.0f days',days);
end
if hours > 0
  str = sprintf('%s%.0f hours',join(str),hours);
end
if mins > 0
  str = sprintf('%s%.0f mins',join(str),mins);
end
if secs > 0
  str = sprintf('%s%.0f secs',join(str),secs);
end

function d=join(str)
  if ~isempty(str)
    d = [str ', '];
  else
    d = [];
  end
end

end