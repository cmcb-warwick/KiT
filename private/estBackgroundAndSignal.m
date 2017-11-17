function [backVal,sigVal,imgHist]=estBackgroundAndSignal(img)
% Estimates value of background as strongest peak and value of signal peak as
% largest peak with value greater than background. Returns NaN for sigVal if
% no signal peak found.
%
% Copyright (c) 2013 Jonathan W. Armond

if ~isvector(img)
  img = img(:);
end

% Find peaks in histogram.
[f,xi]=ksdensity(img,'npoints',100);
imgHist=[xi' f'];

% K-means cluster pixels into background and signal.
[~,c] = kmeans(img,2,'distance','cityblock','emptyaction','singleton');
c = sort(c);
backVal = c(1);
sigVal = c(2);
