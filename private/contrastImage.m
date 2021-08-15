function imgout = contrastImage(imgin,method,range)
%CONTRASTIMAGE  Adjusts the contrast of a 3D image.
%
%   IMGOUT = CONTRASTIMAGE(IMGIN) Normalises elements of matrix IMGIN to
%   fall within the range [0 1] to generate contrast, returning as IMGOUT.
%
%   IMGOUT = CONTRASTIMAGE(IMGIN,METHOD) Specify method of normalisation.
%
%   IMGOUT = CONTRASTIMAGE(IMGIN,METHOD,RANGE) Define range over which the
%   normalisation occurs.
%
%   Options include:
%
%       method: Elements of IMGIN can be normalised using either
%       'percentile' or 'absolute' methods. Default is 'percentile'.
%
%       range:  A two-element vector [a b] where a>=0 and b>a. The
%       'percentile' method requires that b<=1 and a>0, however the
%       'absolute' method has no restriction. Default is [0.1 1].
%
%   Method details:
%
%       'percentile': This method forces all values above b*max(IMGIN) to
%       be 1, all values below a*min(IMGIN) to be 0, and linearly
%       normalises all values in between to be in the range [0 1].
%
%       'absolute': This method forces all values above b to be 1, all
%       values below a to be 0, and linearly normalises all values in
%       between to be in the range [0 1].
%
% Copyright (c) 2018 Chris Smith

% predefine defaults
if nargin<1 || isempty(imgin)
    error('No image given.');
end
if nargin<2 || isempty(method)
    method = 'percentile';
end
if nargin<3 || isempty(range)
    range = [0.1 1];
end

% check input
if ~ismember(method,{'percentile','absolute'})
    error('Unrecognised contrast method: %s',method);
end
if range(1)<0
    error('Error in range: Lower bound must be greater than or equal to zero.')
elseif range(2)<range(1)
    error('Error in range: Upper bound must be greater than lower bound.')
end

% calculate final range
if strcmp(method,'percentile')
    if range(2)>1
        error('Error in range: Percentile method requires that upper bound be no greater than 1.');
    end
	range = stretchlim(imgin,range);
end

% generate output image
imgout = imadjust(imgin,range,[]);

end


