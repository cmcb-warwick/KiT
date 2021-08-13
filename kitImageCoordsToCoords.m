function movieInfo = kitImageCoordsToCoords(dataStruct)
%provide a function to switch between coordinates in our chosen coordinate
%system and image coordinates
%
%Jonathan U Harrison 2019-02-20
%%%%%%%%%%%%%%%%

%get number of time points in movie
nTimepoints = dataStruct.dataProperties.movieSize(4);

%get kinetochore coordinates and amplitude
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),...
    nTimepoints,1);

%if rotated coordinates are to be used ...
if dataStruct.dataProperties.tracksParam.rotate == 1 && ...
        ~isempty(dataStruct.planeFit)

    %get the rotated coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.planeFit(iTime).rotatedCoord;
        if ~isempty(allCoord)
            movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
            movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
            movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
            movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
        end
    end

else %if the original coordinates are to be used
    centerOfMass = zeros(nTimepoints,3);

    %get the original coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.initCoord(iTime).allCoord;
        if ~isempty(allCoord)
            movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
            movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
            movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
            movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;

            %calculate the center of mass in each frame
            centerOfMass(iTime,:) = [mean(movieInfo(iTime).xCoord(:,1)) ...
                mean(movieInfo(iTime).yCoord(:,1)) ...
                mean(movieInfo(iTime).zCoord(:,1))];

            %shift coordinates by center of mass to make the origin in each frame
            %at its center of mass
            movieInfo(iTime).xCoord(:,1) = movieInfo(iTime).xCoord(:,1) - centerOfMass(iTime,1);
            movieInfo(iTime).yCoord(:,1) = movieInfo(iTime).yCoord(:,1) - centerOfMass(iTime,2);
            movieInfo(iTime).zCoord(:,1) = movieInfo(iTime).zCoord(:,1) - centerOfMass(iTime,3);
        end
    end
end