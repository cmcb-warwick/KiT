function [time_strct] = sec2struct(x)
%% function to translate seconds to hh:mm:ss.ms struct
% also creates a string formatted text with the data
% Input: x - scalar value of seconds
% Output: time_strct struct with the fields: 
%                   hour
%                   minutes
%                   seconds
%                   miliseconds
% Date: 5/5/05.
% Created by: Tal Levinger
% email: skyman76@hotmail.com




time_strct.hour = floor(x/60^2);
time_strct.minutes = floor(mod((x/60), 60));
time_strct.seconds = floor(mod(x,60));
time_strct.miliseconds = mod(x,1)*1e3;

time_strct.str = sprintf('%02d:%02d:%02d.%3.0f', time_strct.hour, time_strct.minutes, time_strct.seconds, time_strct.miliseconds);
