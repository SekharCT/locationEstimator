clc;clear;close all;

%% Initialize Variables

% Variables regarding Satellite
% Satellite Velocity 
sysConfig.satelliteVelocity       = 7800; % m/s
% Satellite Distance from Cenre of Earth
sysConfig.distSatelliteCentre     = 7000; % m
% Satellite Location 1
sysConfig.satelliteLoc1           = [38.501889, -121.520728]; % [Latitude, Longitude]
% Satellite Location 2
sysConfig.satelliteLoc2           = [38.596828, -121.531333]; % [Latitude, Longitude]

% Variables regarding UE
% Radius of Earth
sysConfig.earthRadius             = 6400; %m
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

%% Iterative algorithm
while delta < deltaThreshold

    % Get the UE Location
    sysConfig.newUeLoc = getUeLoc(sysConfig);

    % Compute Angle between Norths
    sysConfig.newAngleBetweenNorths = sysConfig.newUeLoc(2) - sysConfig.satelliteLoc2(2);

    % Compute Delta
    delta = sysConfig.newAngleBetweenNorths - sysConfig.oldAngleBetweenNorths;

    % Updaing the Old Angle between the norths
    sysConfig.oldAngleBetweenNorths = sysConfig.newAngleBetweenNorths;

end % End of while loop

%% Function Details
% Name              - getUeLoc
% Details           - Computes User Location using Iterative procedure
% Input Parameters  - sysConfig         - Structure with all the parameters
%                                         needed
% Output Parameters - ueLocation        - User Location
function ueLocation = getUeLoc(sysConfig)
    % Projection Angle of Satellite Velocity and Plane containing User,
    % Satellite and Centre of Earth
    satelliteProjectionAngle = sysConfig.satelliteAngleToNorth + ...
                                sysConfig.oldAngleBetweenNorths + ...
                                sysConfig.ueNorthAdjustedAoA;

    % Achieving Angle of Depression from the satellie
    sysConfig.phi = acos(sysConfig.observedDoppler/(sysConfig.satelliteVelocity*cosd(satelliteProjectionAngle)));

    % Achieving the distance of the Satellite to User
    sysConfig.distSatellite2User = sqrt(sysConfig.earthRadius.^2 + sysConfig.distSatelliteCentre.^2 - ...
                                2*sysConfig.earthRadius.*(sysConfig.distSatelliteCentre.*cos(sysConfig.phi)));

    % Get UE Location using Distance, Angle of Depression and Angle of
    % Arrival
    ueLocation = computeUeCoordinates(sysConfig);
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
% Name              - getSatelliteAngletoNorth
% Details           - Computes Satellite's velocity vector's angle to it's
%                     North
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - satelliteAngletoNorth    - Angle of Arrival
function satelliteAngletoNorth = getSatelliteAngletoNorth(sysConfig)

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