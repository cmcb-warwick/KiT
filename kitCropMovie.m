function [crop,cropSize]=kitCropMovie(movieFileName,zPlane)
%KITCROPMOVIE Presents movie projection for user to crop
%
%  SYNOPSIS crop=kitCropMovie(movieFileName)
%
%  INPUT movieFileName: Filename of movie to crop.
%        zPlane: Show single z-plane (optional).
%
%  OUTPUT crop: Vector of crop coordinates as [xmin,ymin,xmax,ymax]
%
% Copyright (c) 2012 Jonathan W. Armond

if nargin<2
  zPlane=[];
end

try
  [f,rgbImg,zPlanes] = kitMovieProj(movieFileName,zPlane,[],1);
catch me
  errordlg(sprintf('Error opening file %s: %s',movieFileName,me.message));
  crop = []; cropSize = [];
  return
end
fax  = f.CurrentAxes;

% Show image with crop tool.
crop = [];
p = [];
title(fax,'Draw ROI rectangles. Double click to record each.');
finBtn = button(f,'Finish',[2 1 12 2.5],@(hObj,event) close(f),14);
addBtn = button(f,'Add ROI',[16 1 12 2.5],@addROI_cb,14);
fcn = makeConstrainToRectFcn('imrect',get(fax,'XLim'),get(fax,'YLim'));
while ishghandle(f)
  uiwait(f);
  if ~isempty(p) && ishghandle(f)
    rectangle('Parent',fax,'Position',p,'LineWidth',2,'EdgeColor','y');
    crop = [crop; p];
    p = [];
  end
end

crop = round(crop);

% Process ROIs into crop rectangles.
sz = size(rgbImg);
if isempty(crop)
  % If no ROI selected use whole image.
  crop = [1 1 sz(1:2)];
end

for i=1:size(crop,1)
  cropImg = imcrop(rgbImg,crop(i,:));
  sz = size(cropImg);
  cropSize(i,:) = [sz(1:2) zPlanes];
end

function addROI_cb(hObj,eventdata,handles)
  set(addBtn,'Enable','off');
  set(finBtn,'Enable','off');
  h = imrect(fax,'PositionConstraintFcn',fcn);
  p = wait(h);
  delete(h);
  set(addBtn,'Enable','on');
  set(finBtn,'Enable','on');
  uiresume(f);
end

end
