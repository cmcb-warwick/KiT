function k = loc_max3Df(image,mask)
%LOC_MAX3DF finds the local (intensity) maxima in a 3D image fast
%
%SYNOPSIS  locMaxCoord = locmax3d(image,mask)
%
%INPUT     image: any 3D-matrix of at least size [3,3,3]
%          mask (opt): size of the patch in which the local maximum has to
%                      be a maximum. Has to be odd size. Default: [3,3,3]
%
%OUTPUT    locMaxCoord: [x,y,z]-coordinates of the local maxima
%          function ignores border maxima.
%
%REMARKS:  a slower, but more easily readable version of the code can be
%           found commented at the end of the file
%
%created, but not commented by dT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------test input----------------------------

%look for image
if nargin < 1 || isempty(image)
    error('we need at least an image as input argument!')
end

imgSize = size(image);

%image size
if length(imgSize)~=3 || any(imgSize<3)
    error('we need a 3D image of at least size 3x3x3!')
end

%mask
if nargin < 2 || isempty(mask)
    mask = [3 3 3];
else
    %check that we have a valid mask size (3D, not larger than image)

    if length(mask)~=3 || any(mask>imgSize)
        error(['we need a 3D-mask that is smaller than the input image (maskSize: ',num2str(mask),', imgSize: ',num2str(imgSize),')!']);
    end
    if any(floor(mask/2)==mask/2)
        error(['we need a mask of odd size (maskSize: ',num2str(mask),')!']);
    end
end

%-----------end test input-------------------------------


%-----------do doms fast but crazy routine

s = imgSize;
ci=zeros(prod(imgSize),1);
ct = 0;
for p=2:(s(3)-1)
   tmp =  (find((image(2:(s(1)-1),1:(s(2)-2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),3:s(2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),1:(s(2)-2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),3:s(2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),1:(s(2)-2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),3:s(2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),2:(s(2)-1),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),2:(s(2)-1),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),1:(s(2)-2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),1:(s(2)-2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),3:s(2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),3:s(2),p) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),2:(s(2)-1),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),2:(s(2)-1),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),2:(s(2)-1),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),2:(s(2)-1),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),2:(s(2)-1),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),2:(s(2)-1),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),1:(s(2)-2),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),1:(s(2)-2),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),3:s(2),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(2:(s(1)-1),3:s(2),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),1:(s(2)-2),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),3:s(2),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),3:s(2),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),1:(s(2)-2),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),1:(s(2)-2),p-1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),3:s(2),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(1:(s(1)-2),1:(s(2)-2),p+1) < image(2:(s(1)-1),2:(s(2)-1),p))...
        & (image(3:s(1),3:s(2),p-1) < image(2:(s(1)-1),2:(s(2)-1),p)))...
        +(p-2)*(s(2)-2)*(s(1)-2));
    step = length(tmp);
    ci(ct+1:ct+step) = tmp;
    ct = ct + step;
end;

% convert coords of maxima
[l,m,n]=ind2sub(s-2,ci(1:ct));

% add 1, because we ignored border pixels
k=[l+1 m+1 n+1];
% preassign ot
ot = zeros(size(k));

d=floor(mask/2);
ct=1;
for i=1:size(k,1)
    if(all((k(i,:)-d)>0) && all((k(i,:)+d)<=[size(image,1) size(image,2) size(image,3)]))
        patch=image(k(i,1)-d(1):k(i,1)+d(1),k(i,2)-d(2):k(i,2)+d(2),k(i,3)-d(3):k(i,3)+d(3));
        %only the values which are max in mask survive
        if all(image(k(i,1),k(i,2),k(i,3))>=patch)
            ot(ct,:)=k(i,:);
            ct=ct+1;
        end;
    end;
end;
k=ot(1:ct-1,:);
