clc;clear;close all;

%% Initialize Variables

% Variables regarding Satellite
% Satellite Velocity 
sysConfig.satelliteVelocity       = 7800; % m/s
% Satellite Distance from Cenre of Earth
sysConfig.distSatellite2Centre     = 7000; % m
% Satellite Location 1
sysConfig.satelliteLoc1           = [38.501889, -121.520728]; % [Latitude, Longitude]
% Satellite Location 2
sysConfig.satelliteLoc2           = [38.596828, -121.531333]; % [Latitude, Longitude]

% Variables regarding UE
% Radius of Earth
sysConfig.earthRadius             = 6400; % m
% User Location
sysConfig.ueLocation              = [37.389108, -122.143192]; % [Latitude, Longitude] 
% User Azimuth wrt North
sysConfig.ueAzimuthWN             = 30; % Degrees towards East

% AoA Computation
sysConfig.genieBasedAoAComp       = true;

% Doppler Computation
sysConfig.genieBasedDopplerComp   = true;

% Thresholds and Initializations
sysConfig.deltaThreshold          = 0.01;
sysConfig.delta = 1000;
sysConfig.oldAngleBetweenNorths = 0;

%% Calculating Satellite's Angle to North correction for velocity vector
sysConfig.satelliteAngleToNorth = getSatelliteAngletoNorth(sysConfig);

%% Calculating AoA
[sysConfig.observedAoA, sysConfig.ueNorthAdjustedAoA] = computeAoA(sysConfig);

%% Calculating Doppler
sysConfig.observedDoppler = computeDoppler(sysConfig);

%% Fetching UE Location
sysConfig.estimatedUeLocation = getUeLoc(sysConfig);

%% Dispalying Results and Accuracy 
[sysConfig.ueLocation, sysConfig.estimatedUeLocation, ...
    sqrt(sum((convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation) - ...
                convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.estimatedUeLocation)).^2))]

%% Function Details
% Name              - getUeLoc
% Details           - Computes User Location using Iterative procedure
% Input Parameters  - sysConfig         - Structure with all the parameters
%                                         needed
% Output Parameters - ueLocation        - User Location
function ueLocation = getUeLoc(sysConfig)

    % Achieving the distance of the Satellite(S) to User(U)
    [sysConfig.distSatellite2User, sysConfig.angleBetnSCnCU] = getDistUser2Satellite(sysConfig);

    % Angle between velocity vector projected on SCU Plane 
    sysConfig.phi = sysConfig.angleBetnSCnCU + 90;

    % Projection Angle of Satellite Velocity and Plane containing User,
    % Satellite and Centre of Earth
    sysConfig.satelliteProjectionAngle = acosd(sysConfig.observedDoppler/(sysConfig.satelliteVelocity*cosd(sysConfig.phi)));

    % sysConfig.satelliteProjectionAngle = sysConfig.satelliteAngleToNorth + ...
    %                                         sysConfig.oldAngleBetweenNorths + ...
    %                                         sysConfig.ueNorthAdjustedAoA;
    % sysConfig.angleBetnCsuPlaneNorth = sysConfig.oldAngleBetweenNorths + ...
    %                                         sysConfig.ueNorthAdjustedAoA;
    sysConfig.angleBetnCsuPlaneNorth = sysConfig.satelliteProjectionAngle - ...
                                            sysConfig.satelliteAngleToNorth;

    % Iterative algorithm to achieve Angles between Norths;
    while delta < deltaThreshold

        % UE Azimuth Angle wrt North
        sysConfig.currUeFinalAzimuthAngle = sysConfig.angleBetnCsuPlaneNorth - sysConfig.oldAngleBetweenNorths;

        % Get the UE Location
        sysConfig.newUeLoc = computeUeCoordinates(sysConfig);

        % Compute Angle between Norths
        sysConfig.newAngleBetweenNorths = sysConfig.newUeLoc(2) - sysConfig.satelliteLoc2(2);

        % Compute Delta
        delta = sysConfig.newAngleBetweenNorths - sysConfig.oldAngleBetweenNorths;

        % Updaing the Old Angle between the norths
        sysConfig.oldAngleBetweenNorths = sysConfig.newAngleBetweenNorths;

    end % End of while loop
end


%% Function Details
% Name              - computeUeCoordinates
% Details           - Computes UE Location using it's distance, angle of
%                     depression wrt satellite and user's angle of arrival
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - ueLocation     - User's Location Coordinates
function ueLocation = computeUeCoordinates(sysConfig)

end

%% Function Details
% Name              - getDistUser2Satellite
% Details           - Computes Satellite's distance from User
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - distSatellite2User    - Satellite to User Distance
%                   - angleBetnSCnCU        - Angle between line connecting 
%                                             Satellite, Centre of Earth and 
%                                             Satellite, User in Satellite, User 
%                                             and Centre of Earth Plane.
function [distSatellite2User, angleBetnSCnCU] = getDistUser2Satellite(sysConfig)

    % Converting Latitudes and Longitudes to Cartesian Coordinates
    earthCentreCoordinate = [0, 0, 0];
    satelliteCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
    ueCoordinates = convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation);

    % Distance between each points
    % S - Satellite, U - User, C - Centre of Earth
    SC = sqrt(sum((satelliteCoordinates - earthCentreCoordinate).^2));
    UC = sqrt(sum((ueCoordinates - earthCentreCoordinate).^2));

    % We know SC, UC, and angle between UC, SU (angleOfElevation)
    % SU can be solved using the equation solving:
    % SU^2 - 2*SU*UC*cos(angleOfElevation) + UC^2 - SC^2 = 0
    % SU = (2*UC*cosd(angleOfElevation) - sqrt(((2*UC*cosd(angleOfElevation))^2) - ...
    %           4*((UC^2)-(SC^2))))/2
    distSatellite2User = UC*cosd(sysConfig.phi) - sqrt(((UC*cosd(sysConfig.phi)).^2) - ...
                            ((UC.^2) - (SC.^2)));
end

%% Function Details
% Name              - getSatelliteAngletoNorth
% Details           - Computes Satellite's velocity vector's angle to it's
%                     North
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - satelliteAngletoNorth    - Angle of Arrival
function satelliteAngletoNorth = getSatelliteAngletoNorth(sysConfig)
    % Consider the Locations to fetch the angle
    velocityVectorStart = sysConfig.satelliteLoc1;
    velocityVectorEnd = sysConfig.satelliteLoc2;
    northAxisStart = sysConfig.satelliteLoc1;
    northAxisEnd = [sysConfig.satelliteLoc2(1), sysConfig.satelliteLoc1(2)];

    % Converting to Cartesian Vectors
    velocityVector = convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorEnd) - ...
                        convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorStart);
    northAxisVector = convLatLong2Cartesian(sysConfig.distSatellite2Centre, northAxisEnd) - ...
                        convLatLong2Cartesian(sysConfig.distSatellite2Centre, northAxisStart);

    % Angle between the vectors
    satelliteAngletoNorth = acosd(sum(velocityVector.*northAxisVector)/(sqrt((sum(velocityVector.^2))*(sum(northAxisVector.^2)))))
end

%% Function Details
% Name              - computeDoppler
% Details           - Computes Doppler
% Input Parameters  - sysConfig         - Structure with all the parameters
%                                         needed
% Output Parameters - observedDoppler   - Observed Doppler
function observedDoppler = computeDoppler(sysConfig)

end

%% Function Details
% Name              - computeAoA
% Details           - Computes Angle of Arrival
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - observedAoA    - Angle of Arrival
function observedAoA = computeAoA(sysConfig)
    % Converting Latitudes and Longitudes to Cartesian Coordinates
    earthCentreCoordinate = [0, 0, 0];
    satelliteCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
    ueCoordinates = convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation);

    % Distance between each points
    % S - Satellite, U - User, C - Centre of Earth
    SC = sqrt(sum((satelliteCoordinates - earthCentreCoordinate).^2));
    UC = sqrt(sum((ueCoordinates - earthCentreCoordinate).^2));
    SU = sqrt(sum((satelliteCoordinates - ueCoordinates).^2));

    % Angle of Arrival
    observedAoA = acosd(((UC.^2) + (SU.^2) - (SC.^2))./(2*(UC.*SU)));
end


%% Function Details
% Name              - convLatLong2Cartesian
% Details           - Converts Latitude, Longitude domain to cartesian
%                   - domain
% Input Parameters  - radius                - Radius from the centre
%                   - latLongPair           - Latitude, Longitude Value
% Output Parameters - cartesianCoordinates  - [x,y,z] Coordinates
function cartesianCoordinates = convLatLong2Cartesian(radius, latLongPair)

    % Filling the coordinates
    x = radius*cosd(latLongPair(1))*cosd(latLongPair(2));
    y = radius*cosd(latLongPair(1))*sind(latLongPair(2));
    z = radius*sind(latLongPair(1));

    % Putting them together
    cartesianCoordinates = [x, y, z];
end


%% Function Details
% Name              - convCartesian2LatLong
% Details           - Converts cartesian doamin to Latitude, Longitude 
%                   - domain
% Input Parameters  - radius                - Radius from the centre
%                   - cartesianCoordinates  - [x,y,z] Coordinates
% Output Parameters - latLongPair           - Latitude, Longitude Value
function latLongPair = convCartesian2LatLong(radius, cartesianCoordinates)

    % Filling the coordinates
    latitude = asind(cartesianCoordinates(3)/radius);
    longitude = acosd(cartesianCoordinates(1)/radius/cosd(latitude));

    % Putting them together
    latLongPair = [latitude, longitude];
end