function stack=kitReadImageStack(imageReader,metadata,t,c,crop,normalize)
% KITREADIMAGESTACK Read a single image frame stack from a movie file
%
%    STACK = KITREADIMAGESTACK(IMAGEREADER,METADATA,T,C,CROP,NORMALIZE) Read
%    a single image frame stack (i.e. all z-planes) at time T in channel C from
%    IMAGEREADER described by METADATA.
%
%    CROP Optional, vector of [XMIN,YMIN,WIDTH,HEIGHT] for cropping stack.
%
%    NORMALIZE Optional, 0, 1 or -1. Normalize by maximum pixel value. Defaults
%    to 1. If -1, no normalization is performed and image is returned in
%    original datatype, otherwise it is converted to double.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<5
  crop = [];
end

if nargin<6
  normalize = 0;
end

if normalize == -1
  dataType = metadata.dataType;
else
  dataType = 'double';
end

stackSize = kitComputeStackSize(crop,metadata.frameSize);

stack = zeros(stackSize, dataType);
for z = 1:metadata.frameSize(3)
  stack(:,:,z) = kitReadImagePlane(imageReader, metadata, t, c, z, crop, normalize);
end

