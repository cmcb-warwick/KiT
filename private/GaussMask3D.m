function gauss=GaussMask3D(sigma,fSze,cent,cnorm,absCenter,noInt,separated)
% GaussMask3D	create a gaussian 3D mask
%
%    SYNOPSIS gauss =GaussMask3D(sigma,fSze,cent,cnorm,absCenter,separated);
%
%    INPUT: sigma  of gauss mask [sigmaX sigmaY sigmaZ] or [sigma].
%                  If scalar, sigma will be the same for all dimensions
%           fSze   size of the gauss mask [sizeX sizeY sizeZ] or [size]
%                  If scalar, size will be the same for all dimensions.
%                  (odd size required for symmetric mask!)
%           cent   (optional)3D vector with center position [0 0 0] is
%                  center of fSze (default=[0 0 0])
%                  if absCenter = 1, center position is measured from
%                  [0,0,0] in matrix coordinates (which is outside of the
%                  array by half a pixel!)
%           cnorm  (optional) select normalization method:
%                  =0 (default) no normalization - max of gauss will be 1
%                  =1 norm so that integral of all pixels will be exactly 1
%                  =2 norm so that integral of an infinite Gauss would be 1
%           absCenter (optional) select the type of center vector
%                  =0 (default). Zero is at center of array
%                  =1 Zero is at [0,0,0] in matrix coordinates, which is
%                  outside of the array by half a pixel!
%                      (center of top left pixel of 2-D matrix = [1,1])
%           noInt  (optional) normally, GaussMask3D would calculate
%                  integral over a voxel. However, due to numerical
%                  problems, this integral will become zero after 8 sigma.
%                  =0 (default) Intensities are the integral over voxels
%                  =1 Intensities are the value of the Gauss at the voxel
%                     center
%           separated (optional; [{0}/1]) if 1, gauss will be a cell array
%                  with 3 vectors, gx, gy, gz that can be used to
%                  sequentially filter the image in order to speed up
%                  calculations
%
%    OUTPUT: gauss   3D gaussian intensity distribution with voxel values
%                    equal to the integral of the gauss over the area of
%                    the voxel
%
% c: 8/05/01 dT
% corrected integral by jonas

% test input
ls = length(sigma);
switch ls
    case 3
        % all is good
    case 1
        sigma = [sigma,sigma,sigma];
    otherwise
        error('sigma has to be either a 1-by-3 vector or a scalar!')
end

lf = length(fSze);
switch lf
    case 3
        % all is good
    case 1
        fSze = [fSze,fSze,fSze];
    otherwise
        error('fSze has to be either a 1-by-3 vector or a scalar!')
end
if nargin < 3 || isempty(cent)
    cent=[0 0 0];
end;
if nargin<4 || isempty(cnorm)
    cnorm = false;
end;
if nargin < 5 || isempty(absCenter)
    absCenter = false;
end
if nargin < 6 || isempty(noInt)
    noInt = false;
end
if nargin < 7 || isempty(separated)
    separated = false;
end

% transform absolute center into relative center
if absCenter
    cent = cent-(fSze+1)/2;
end


% to get the accurate pixel intensities, use the cumsum of the gaussian
% (taken from normcdf.m). As the origin is in the center of the center
% pixel, the intensity of pixel +2 is the integral from 1.5:2.5, which,
% fortunately, has the spacing 1.

gauss=zeros(fSze);
x=((-fSze(1)/2:fSze(1)/2)-cent(1))./sigma(1);
y=((-fSze(2)/2:fSze(2)/2)-cent(2))./sigma(2);
z=((-fSze(3)/2:fSze(3)/2)-cent(3))./sigma(3);
if noInt
    % calculate the Gaussian
    ex = 1/sqrt(pi*2) * exp(-(x(1:end-1)+0.5/sigma(1)).^2./2)/sigma(1);
    ey = 1/sqrt(pi*2) * exp(-(y(1:end-1)+0.5/sigma(2)).^2./2)/sigma(2);
    ez = 1/sqrt(pi*2) * exp(-(z(1:end-1)+0.5/sigma(3)).^2./2)/sigma(3);
else
    % calculate the integral
    ex = diff(0.5 * erfc(-x./sqrt(2)));
    ey = diff(0.5 * erfc(-y./sqrt(2)));
    ez = diff(0.5 * erfc(-z./sqrt(2)));
end

if separated
    % return cell array of vectors

    % norm Gauss
    switch cnorm
        case 0 % maximum of Gauss has to be 1
            ex = ex*((2*pi)^0.5*prod(sigma));
            ey = ey*((2*pi)^0.5*prod(sigma));
            ez = ez*((2*pi)^0.5*prod(sigma));

        case 1
            ex = ex/sum(ex);
            ey = ey/sum(ey);
            ez = ez/sum(ez);

        case 2 % the whole erfc thing is already normed for infinite Gauss
            % so nothing to do here.
            % gauss(:) = gauss;

        otherwise
            % no change to gauss

    end

    gauss = {ex',ey,permute(ez,[3,1,2])};

else
    % construct the 3D matrix (nice work by Dom!)
    exy=ex'*ey;
    gauss(:)=exy(:)*ez;

    % norm Gauss
    switch cnorm
        case 0 % maximum of Gauss has to be 1
            gauss = gauss*((2*pi)^1.5*prod(sigma));

        case 1
            gauss = gauss/sum(gauss(:));
        case 2 % the whole erfc thing is already normed for infinite Gauss
            % so nothing to do here.
            % gauss(:) = gauss;

        otherwise
            % no change to gauss

    end
end
