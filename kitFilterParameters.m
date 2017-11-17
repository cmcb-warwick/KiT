function [sizeXY, sizeZ]=kitFilterParameters(wavelength, na, pixelSize)
%KITFILTERPARAMETERS Define PSF filter parameters
%
% SYNOPSIS: job=kitFilterParameters(job)
%
% INPUT wavelength: wavelength in microns (optional, default 0.525)
%
%       na: numerical aperture (optional, default 1.4)
%
%       pixelSize : [sizeX, sizeZ]. If size is given, the
%       psf-width is returned in pixel (optional, default [1,1])
%
%
% OUTPUT sizeXY, sizeZ: Size of psf in microns or pixels.
%
% Copyright (c) 2007 Jonas Dorn
% Copyright (c) 2012 Jonathan W. Armond

if nargin<1
    wavelength = 0.525;
end
if nargin<2
    na = 1.4;
end
if nargin<3
    pixelSize = [1,1];
end

% Refractive index of oil and air.
n = 1.51;

% multiplication factor for Gauss fitted to Bessel (see Dom's
% first JM paper)
gob = [0.21, 0.66];

% Calculate psf-size
sizeXY = (gob(1)*wavelength/na)/pixelSize(1);
sizeZ  = (gob(2)*n*wavelength/na^2)/pixelSize(2);