function [spots,spotsPix,img]=kitFrameSpots(filename,channel)

job = kitDefaultOptions();
job.movie = filename;
[crop, cropSize] = kitCropMovie(job.movie);
job.crop = crop;
job.cropSize = cropSize;

[job.metadata,reader] = kitOpenMovie(filename);

%hack metadata
job.metadata.pixelSize = [0.0694 0.0694 0.1];
job.metadata.wavelength(1) = 0.647;

job.options.debug.showMmfCands = 1;
job.options.debug.showMmfClusters = 0;
job.options.debug.showMmfFinal = -2;
job.options.debug.showCentroidFinal = -2;
job.options.debug.mmfVerbose = 1;
%job.options.clusterSeparation = 1;
job.options.alphaF = 0;
job.options.alphaA = 0.5;

%dataStruct = kitMakeMakiDatastruct(job, channel);
%dataStruct = kitFindCoordsMaki(job, reader, dataStruct, channel);
%job.dataStruct{channel} = dataStruct;
job.dataStruct{channel} = kitMakeMakiDatastruct(job, channel);
job = kitMixtureModel(job,reader,channel);

spots = job.dataStruct{channel}.initCoord(1).allCoord(:,1:3);
spotsPix = job.dataStruct{channel}.initCoord(1).allCoordPix(:,1:3);

img = kitReadImageStack(reader,job.metadata,1,channel,job.crop,1);
