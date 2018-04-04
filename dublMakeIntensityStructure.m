function intraStruct = dublMakeIntensityStructure()
% DUBLMAKEINTENSITYSTRUCTURE Produces an empty structure into which
% kinetochore intensity measurements are compiled using
% dublIntensityMeasurements.
%
%    DUBLMAKEINTENSITYSTRUCTURE() Produces a structure allowing for
%    kinetochore measurements to be compiled in order to draw
%    population-scale analyses. No input is required.
%
%
% Copyright (c) 2018 C. A. Smith

intraStruct.kitVersion = kitVersion;
intraStruct.label = [];

% make substructures for all pair-derived measurements
subint = struct('inner',[],'outer',[]);
intensity = struct('mean',subint,'max',subint,'bg',subint);

intraStruct.intensity = intensity;
    
end
    
    
    
    