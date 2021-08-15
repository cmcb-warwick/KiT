function interStruct = kitMakeInterStructure
% KITMAKEINTERSTRUCTURE Produces an empty structure into which
% inter-kinetochore measurements are compiled using kitInterMeasurements.
%
%    KITMAKEINTERSTRUCTURE() Produces a structure allowing for inter-
%    kinetochore measurements to be compiled in order to draw
%    population-scale analyses. No input is required.
%
%
% Copyright (c) 2018 C. A. Smith


interStruct.kitVersion = kitVersion;
interStruct.label = [];

% make substructures for all pair-derived measurements
direction = struct('P',[],'AP',[],'N',[],'S',[]);
sisSep = struct('x',[],'y',[],'z',[],...
               'twoD',[],'threeD',[]);
intensity = struct('mean',[],'max',[],'bg',[]);

coords = struct('x',[],'y',[],'z',[]);
% those only appropriate for plate coordinate system
twist  = struct('y',[],'z',[],'threeD',[]);
sisterCentreSpeed = [];
plateThickness = [];

% produce generic substructure for microscope coordinate system
microscope.coords = coords;
microscope.sisSep = sisSep;
% produce same again for plate coordinate system
plate = microscope;
plate.twist = twist;
plate.sisterCentreSpeed = sisterCentreSpeed;
plate.plateThickness    = plateThickness;
  
interStruct.direction = direction;
interStruct.microscope = microscope;
interStruct.plate = plate;
interStruct.intensity = intensity;


end
    
    
    
    