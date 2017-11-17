function movie=kitReadWholeMovie(imageReader,metadata,c,crop,normalizePlanes,normalize)
% KITREADIMAGESTACK Read a whole movie from a movie file
%
%    MOVIE = KITREADIMAGESTACK(IMAGEREADER,METADATA,C,CROP,NORMALIZE) Read a
%    whole movie in channel C from IMAGEREADER described by METADATA.
%
%    CROP Optional, vector of [XMIN,YMIN,WIDTH,HEIGHT] for cropping stack.
%
%    NORMALIZEPLANES Optional, 0, 1 or -1. Normalize by maximum pixel value. Defaults
%    to 0. If -1, no normalization is performed and image is returned in
%    original datatype, otherwise it is converted to double.
%
%    NORMALIZE Optional, 0 or 1. Normalize by maximum pixel entire movie. Defaults to 0.
%
%    Alternatively, IMAGEREADER can be a string filename as a shortcut.
%
% Copyright (c) 2013 Jonathan W. Armond
if nargin<3
  c = 1;
end

if nargin<4
  crop = [];
end

if nargin<5
  normalizePlanes = 0;
end
if nargin<6
  normalize = 0;
end

if ischar(imageReader)
  % imageReader is a filename.
  [metadata,rdr] = kitOpenMovie(imageReader);
  movie = kitReadWholeMovie(rdr,metadata,c,crop,normalizePlanes,normalize);
  rdr.close();
  return
end

if normalizePlanes == -1
  dataType = metadata.dataType;
else
  dataType = 'double';
end

stackSize = kitComputeStackSize(crop,metadata.frameSize);

movie = zeros([stackSize, metadata.nFrames], dataType);
for t = 1:metadata.nFrames
  movie(:,:,:,t) = kitReadImageStack(imageReader, metadata, t, c, crop, normalizePlanes);
end

if normalize>0
  movie = movie/max(movie(:));
end
