function [icov,centroid]=imageCov(img)
% IMAGECOV Calculate image covariance matrix and centroid of image.
%
% created: jarmond (08/2012)

% Y refers to ROWS, X refers to COLS in IMG, hence this ordering.
[sy,sx,sz] = size(img);
[y,x,z] = ndgrid(1:sy, 1:sx, 1:sz);

% Raw moments.
M_000 = sum(img(:));
i = x.*img;
M_100 = sum(i(:));
i = y.*img;
M_010 = sum(i(:));
if sz>1
    i = z.*img;
    M_001 = sum(i(:));
end

i = (x.^2).*img;
M_200 = sum(i(:));
i = (y.^2).*img;
M_020 = sum(i(:));
if sz>1
    i = (z.^2).*img;
    M_002 = sum(i(:));
end

i = x.*y.*img;
M_110 = sum(i(:));
if sz>1
    i = x.*z.*img;
    M_101 = sum(i(:));
    i = y.*z.*img;
    M_011 = sum(i(:));
end

% Central moments.

if sz>1
    centroid = [M_100 M_010 M_001]./M_000;
else
    centroid = [M_100 M_010]./M_000;
end

u_200 = M_200 - centroid(1)*M_100;
u_020 = M_020 - centroid(2)*M_010;
if sz>1
    u_002 = M_002 - centroid(3)*M_001;
end

u_110 = M_110 - centroid(1)*M_010;
if sz>1
    u_101 = M_101 - centroid(1)*M_001;
    u_011 = M_011 - centroid(2)*M_001;
end

if sz>1
    icov = [u_200 u_110 u_101;
            u_110 u_020 u_011;
            u_101 u_011 u_002]./M_000;
else
    icov = [u_200 u_110;
            u_110 u_020]./M_000;    
end