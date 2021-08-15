function [imgCrpd,coords] = indivSpot(img,coords,opts)

% get pixel resolution
imageSize = size(img);
imgw = opts.imgHalfWidth;
bgcol = opts.bgcol;

% predefine cropped image
imgCrpd = ones(2*imgw+1)*bgcol;

% calculate centre pixel
centrePxl = round(coords);

% max project over three z-slices around point
img = max(img(:,:,max(1,centrePxl(3)-2):min(centrePxl(3)+2,opts.imageSize(3))), [], 3);

% produce cropped image around track centre
xReg = [max(centrePxl(1)-imgw+1,1) min(centrePxl(1)+imgw+1,imageSize(2))];
yReg = [max(centrePxl(2)-imgw+1,1) min(centrePxl(2)+imgw+1,imageSize(1))];
imgCrpd(1:diff(yReg)+1,1:diff(xReg)+1) = img(yReg(1):yReg(2),xReg(1):xReg(2));

% define contrast stretch and apply
irange=stretchlim(imgCrpd,[0.1 1]);
imgCrpd = imadjust(imgCrpd, irange, []);

% correct coordinates to the cropped region
coords(:,1:2) = coords(:,1:2) - [xReg(1) yReg(1)];

end

