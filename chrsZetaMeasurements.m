function compiledZeta = chrsZetaMeasurements(images,varargin)
% DUBLINTRAMEASUREMENTS Produces a structure of population-level
% intra-kinetochore measurements over multiple experiments.
%
%    CHRSZETAMEASUREMENTS(IMAGES,...) Provides all coordinates and inter-
%    probe measurements for all spots across all experiments used for
%    chromatic shift calculation.  The resulting structure provides the
%    tools for deriving population-level analyses for an experiment.
%    Options are available.
%
%    Options, defaults in {}:-
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. The channels between
%       which to make intra-kinetochore measurements. The direction of
%       measurements will be defined by the channel orientation in the
%       neighbourSpots section of options.
%
%    prevMeas: {[]} or a structure previously generated. A structure of
%       results from previous experiments to allow new experiment data to
%       be appended.
%
%
% Copyright (c) 2017 C. A. Smith

% default options
opts.channels = [1 2];
opts.prevMeas = [];
% user options
opts = processOptions(opts,varargin{:});

%% Pre-processing input structure

%check structure of movies
if ~iscell(images{1})
    images = {images};
    kitLog('Movie structure provided implies only one experiment. Assuming only one experiment.');
end
%find number of movies
nExpts = length(images);
% find channel vector
chanVect = opts.channels;

%% Preprocessing output structure

if isempty(opts.prevMeas)
    
    % make new intra-measurements structure if no previous measurements
    % provided
    allData = chrsMakeZetaStructure;
    allData.channelVector = opts.channels;
    allData = struct2strForm(allData);
    
else
    % get all old data
    allData = struct2strForm(opts.prevMeas);
    
end

% predesignate error arrays
noDS = [];
noSpot = [];

%% Compiling measurements

for iExpt = 1:nExpts
    
    % get movie and sister list
    theseImages = images{iExpt};
    
    for iMov = 1:length(theseImages)
      
      % get the movie index
      movNum = theseImages{iMov}.index;
      % check whether there is data in this movie
      if ~isfield(theseImages{iMov},'dataStruct')
        noDS = [noDS; iExpt movNum];
        continue
      end
      
      % get dataStructs
      dStrue = theseImages{iMov}.dataStruct{chanVect(1)};
      dSshifted = theseImages{iMov}.dataStruct{chanVect(2)};
      
      % check whether the movie failed
      if ~isfield(dStrue,'failed') || dStrue.failed || ~isfield(dSshifted,'failed') || dSshifted.failed
        noSpot = [noSpot; iExpt movNum];
        continue
      end
        
      % start counter for storing data
      c=1;

      % get dataStructs
      dStrue = theseImages{iMov}.dataStruct{chanVect(1)};
      dSshifted = theseImages{iMov}.dataStruct{chanVect(2)};
      % get raw initCoords
      if isfield(dStrue,'rawInitCoord')
        iCtrue = dStrue.rawInitCoord(1);
        iCshifted = dSshifted.rawInitCoord(1);
      else
        iCtrue = dStrue.initCoord(1);
        iCshifted = dSshifted.initCoord(1);
      end

      %% Kinetochore positions

      % predefine various variables
      micrCoordsTrue = []; micrCoordsShifted = [];
      micrCoord_x = []; micrCoord_y = []; micrCoord_z = [];

      % get microscope coordinates of each spot
      micrCoordsTrue = [micrCoordsTrue; iCtrue.allCoord(:,1:3)];
      micrCoordsShifted = [micrCoordsShifted; iCshifted.allCoord(:,1:3)];
      micrCoord_x = [micrCoord_x; [iCtrue.allCoord(:,1) iCshifted.allCoord(:,1)]];
      micrCoord_y = [micrCoord_y; [iCtrue.allCoord(:,2) iCshifted.allCoord(:,2)]];
      micrCoord_z = [micrCoord_z; [iCtrue.allCoord(:,3) iCshifted.allCoord(:,3)]];

      % put data into string format
      newData(c,:) = {'coords.x',micrCoord_x}; c=c+1;
      newData(c,:) = {'coords.y',micrCoord_y}; c=c+1;
      newData(c,:) = {'coords.z',micrCoord_z}; c=c+1;

      %% Raw zeta measurements
      
      micrData = soloMeasurements(micrCoordsTrue,micrCoordsShifted);
      % put data into string format
      newData(c,:) = {'raw.zeta.x',micrData.zeta_x}; c=c+1;
      newData(c,:) = {'raw.zeta.y',micrData.zeta_y}; c=c+1;
      newData(c,:) = {'raw.zeta.z',micrData.zeta_z}; c=c+1;
      newData(c,:) = {'raw.zeta.twoD',micrData.zeta_2D}; c=c+1;
      newData(c,:) = {'raw.zeta.threeD',micrData.zeta_3D}; c=c+1;

      %% Filtered zeta measurements
      
      % get filtered initCoords
      iCtrue = dStrue.initCoord(1);
      iCshifted = dSshifted.initCoord(1);
      
      % redefine various variables
      micrCoordsTrue = []; micrCoordsShifted = [];
      micrCoord_x = []; micrCoord_y = []; micrCoord_z = [];

      % get microscope coordinates of each spot
      micrCoordsTrue = [micrCoordsTrue; iCtrue.allCoord(:,1:3)];
      micrCoordsShifted = [micrCoordsShifted; iCshifted.allCoord(:,1:3)];
      micrCoord_x = [micrCoord_x; [iCtrue.allCoord(:,1) iCshifted.allCoord(:,1)]];
      micrCoord_y = [micrCoord_y; [iCtrue.allCoord(:,2) iCshifted.allCoord(:,2)]];
      micrCoord_z = [micrCoord_z; [iCtrue.allCoord(:,3) iCshifted.allCoord(:,3)]];
      
      micrData = soloMeasurements(micrCoordsTrue,micrCoordsShifted);
      % put all data into string format
      newData(c,:) = {'filtered.zeta.x',micrData.zeta_x}; c=c+1;
      newData(c,:) = {'filtered.zeta.y',micrData.zeta_y}; c=c+1;
      newData(c,:) = {'filtered.zeta.z',micrData.zeta_z}; c=c+1;
      newData(c,:) = {'filtered.zeta.twoD',micrData.zeta_2D}; c=c+1;
      newData(c,:) = {'filtered.zeta.threeD',micrData.zeta_3D}; c=c+1;

      % compile new data with original
      allData = combineStrForms(allData,newData);  
            
    end % images     
end % expts

%% Save results to structure

compiledZeta = strForm2struct(allData);

%% Output any error information

if ~isempty(noDS)
  fprintf('\nThe following cells failed during spot detection:\n');
  for iCell = 1:size(noDS,1)
    fprintf('    Exp %i, Mov %i\n',noDS(iCell,1),noDS(iCell,2));
  end
end
if ~isempty(noSpot)
  fprintf('\nThe following cells found no spots:\n');
  for iCell = 1:size(noSpot,1)
    fprintf('    Exp %i, Mov %i\n',noSpot(iCell,1),noSpot(iCell,2));
  end
end


function measurements = soloMeasurements(coordsTrue,coordsShifted)
    
    % coordinate-specific delta
    zeta_xyz = coordsShifted(:,1:3) - coordsTrue(:,1:3);
    measurements.zeta_x = zeta_xyz(:,1);
    measurements.zeta_y = zeta_xyz(:,2);
    measurements.zeta_z = zeta_xyz(:,3);
    
    % 3D delta
    zeta_3D = sqrt(sum(zeta_xyz.^2,2));
    measurements.zeta_3D = zeta_3D;
    
    % 2D delta
    zeta_2D = sqrt(sum(zeta_xyz(:,1:2).^2,2));
    measurements.zeta_2D = zeta_2D;

end %soloMeasurements subfunction

end