function [out,M]=Gauss2D(x,sigma);
% Gauss2D	apply a 2 dimensional gauss filter
%
%    out = Gauss2D(x,sigma);
%
%    INPUT: x      image
%           sigma  of gauss filter
%
%    OUTPUT: out   filtered image
%            M     gaussian mask
%
%This file is part of u-track.
%
%    u-track is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%    u-track is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with u-track.  If not, see <http://www.gnu.org/licenses/>.
%
%
%
% Copyright: Ponti: 07/02

R = ceil(3*sigma);   % cutoff radius of the gaussian kernel
M = zeros(2*R+1); % KJ
for i = -R:R,
   for j = -R:R,
      M(i+R+1,j+R+1) = exp(-(i*i+j*j)/2/sigma/sigma);
   end
end
M = M/sum(M(:));   % normalize the gaussian mask so that the sum is
                   % equal to 1
                   
% more correct version - and probably a bit faster
% M = GaussMask2D(sigma,2*R+1,[],1);

% Convolute matrices
out=filter2(M,x);

