function imageCoords = kitCoordsToImageCoords(job,channel,coords,frameNum)
%provide a function to switch between image coordinates
% and coordinatesin our chosen coordinate system
%
% job - struct containing job info including plane fit
%channel - integer for which channel is of interest
%coords - nx3 (or nx6) matrix of coordinates in COM coords
%frameNum - integer from 1 to number of frames in movie
%
%Assume that coords are in rotated frame of reference and scaled by pixel
%size
%Jonathan U Harrison 2019-02-20
%%%%%%%%%%%%%%%%

if ~(isfield(job,'dataStruct') && ...
        isfield(job.dataStruct{channel},'planeFit') && ...
        isfield(job,'metadata') && isfield(job,'options'))
    error('Missing information in job struct');
end

ndims = 2 + job.metadata.is3D;
pixelSize = job.metadata.pixelSize;
chrShift = job.options.chrShift.result{job.options.coordSystemChannel,...
    channel};

if isempty(job.dataStruct{channel}.planeFit(frameNum).planeVectors)
    error('Plane fit needed to convert coordinate systems but failed for this job');
end

%undo previous rotation of the coordinates
unrot = (job.dataStruct{channel}.planeFit(frameNum).planeVectors*...
   coords(:,1:3)')';
%apply chromatic shift and stretch by pixel size
scaled = (unrot + repmat(chrShift(1:ndims), ...
 job.dataStruct{channel}.initCoord(frameNum).nSpots,1)) .* ...
 repmat(1./pixelSize,job.dataStruct{channel}.initCoord(frameNum).nSpots,1);
%shift by centre of mass
centerOfMass = mean(job.dataStruct{channel}.initCoord(frameNum).allCoordPix);
imageCoords = scaled + repmat(centerOfMass(1:3),size(scaled,1),1);

end
