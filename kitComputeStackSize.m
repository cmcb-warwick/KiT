function stackSize=kitComputeStackSize(crop,frameSize)
% KITCOMPUTESTACKSIZE Compute pixel size after cropping

if isempty(crop)
    stackSize = frameSize(1:3);
else
  w = crop(4);
  h = crop(3);
  x1 = crop(2);
  y1 = crop(1);
  x2 = min(round(x1 + w),frameSize(1));
  y2 = min(round(y1 + h),frameSize(2));
  x1 = max(round(x1),1);
  y1 = max(round(y1),1);
  sx = x2-x1+1;
  sy = y2-y1+1;
  stackSize = [sx,sy,frameSize(3)];
end
