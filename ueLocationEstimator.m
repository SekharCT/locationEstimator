clc;clear;close all;

%% Initialize Variables

% Variables regarding Satellite
% Satellite Velocity 
sysConfig.satelliteVelocity       = 7800; % m/s
% Satellite Distance from Cenre of Earth
sysConfig.distSatellite2Centre    = 7000; % m
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

% Movement Sign of SNR
% +1 -> Satellite Going Away or SNR historical trend for this satellite is dereasing.
% -1 -> Satellite Coming Towards or SNR historical trend for this satellite is increasing.
sysConfig.movementSign            = +1;

% AoA Computation
sysConfig.genieBasedAoAComp       = true;

% Doppler Computation
sysConfig.genieBasedDopplerComp   = true;

% Thresholds and Initializations
sysConfig.deltaThreshold          = 0.01;
sysConfig.oldAngleBetweenNorths   = 0;

%% Calculating Satellite's Angle to North correction for velocity vector
sysConfig.satelliteAngleToNorth = getSatelliteAngletoNorth(sysConfig);

%% Calculating AoA
sysConfig.observedAoA = computeAoA(sysConfig);

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
    [sysConfig.distSatellite2User, sysConfig.angleBetnSCnSU] = getDistUser2Satellite(sysConfig);

    % Angle between velocity vector projected on SCU Plane 
    sysConfig.phi = 90 - sysConfig.movementSign*sysConfig.angleBetnSCnSU;

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

    % Iterative algorithm to achieve Angles between Norths
    delta = 10000;
    while delta > sysConfig.deltaThreshold

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

    ueLocation = sysConfig.newUeLoc;
end


%% Function Details
% Name              - computeUeCoordinates
% Details           - Computes UE Location using it's distance, angle of
%                     depression wrt satellite and user's angle of arrival
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - ueLocation     - User's Location Coordinates
function ueLocation = computeUeCoordinates(sysConfig)
    
    % Sign of AoA specified tthe direction of CUS lane with Velocity Vector
    % If +ve -> UE is considered to be on the right of Satellite and vice versa
    if sysConfig.observedAoA > 0
        leftRightSign = +1;
    else
        leftRightSign = -1;
    end
    
    % Direction of velocity projected on CUS Plane
    velocityVectorStart = sysConfig.satelliteLoc1;
    velocityVectorEnd = sysConfig.satelliteLoc2;
    velocityVector = convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorEnd) - ...
                        convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorStart);
    tangentialVelocityVecPlane = convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorStart);
    velocityPerpVectorInCusPlane = [tangentialVelocityVecPlane(2)*velocityVector(3) - ...
                                        velocityVector(2)*tangentialVelocityVecPlane(2), ...
                                        tangentialVelocityVecPlane(3)*velocityVector(1) - ...
                                        velocityVector(3)*tangentialVelocityVecPlane(1), ...
                                        tangentialVelocityVecPlane(1)*velocityVector(2) - ...
                                        velocityVector(1)*tangentialVelocityVecPlane(2)];
    velocityVectorProjectionInCusPlane = velocityVector*cosd(leftRightSign*sysConfig.satelliteProjectionAngle) + ...
                    sind(leftRightSign*sysConfig.satelliteProjectionAngle)*...
                    velocityPerpVectorInCusPlane./sqrt(sum(velocityPerpVectorInCusPlane.^2));

    % Direction on projected velocity in the path direction
    satelliteCentreVector = tangentialVelocityVecPlane;
    orthoCusPlaneVector = [satelliteCentreVector(2)*velocityVectorProjectionInCusPlane(3) - ...
                            velocityVectorProjectionInCusPlane(2)*satelliteCentreVector(2), ...
                            satelliteCentreVector(3)*velocityVectorProjectionInCusPlane(1) - ...
                            velocityVectorProjectionInCusPlane(3)*satelliteCentreVector(1), ...
                            satelliteCentreVector(1)*velocityVectorProjectionInCusPlane(2) - ...
                            velocityVectorProjectionInCusPlane(1)*satelliteCentreVector(2)];
    pathDirection = velocityVectorProjectionInCusPlane*cosd(sysConfig.phi) + ...
                        sind(sysConfig.phi)*orthoCusPlaneVector./sqrt(sum(orthoCusPlaneVector.^2));
    pathDirectionUnitNorm = pathDirection./sqrt(sum(pathDirection.^2));

    % Finding the UE Location in the path direction from Satellite Location
    ueLocationCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorStart) + ...
                               pathDirectionUnitNorm*sysConfig.distSatellite2User;

    % Converting UE Location to Latitude and Longitude Co-ordinates
    ueLocation = convCartesian2LatLong(sysConfig.earthRadius, ueLocationCoordinates);
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
function [distSatellite2User, angleBetnSCnSU] = getDistUser2Satellite(sysConfig)

    % Converting Latitudes and Longitudes to Cartesian Coordinates
    earthCentreCoordinate = [0, 0, 0];
    satelliteCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
    ueCoordinates = convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation);

    % Distance between each points
    % S - Satellite, U - User, C - Centre of Earth
    SC = sqrt(sum((satelliteCoordinates - earthCentreCoordinate).^2));
    UC = sqrt(sum((ueCoordinates - earthCentreCoordinate).^2));

    % Angle of Elevation
    angleOfElevation = abs(sysConfig.observedAoA) + 90;

    % We know SC, UC, and angle between UC, SU (angleOfElevation)
    % SU can be solved using the equation solving:
    % SU^2 - 2*SU*UC*cos(angleOfElevation) + UC^2 - SC^2 = 0
    % SU = (2*UC*cosd(angleOfElevation) - sqrt(((2*UC*cosd(angleOfElevation))^2) - ...
    %           4*((UC^2)-(SC^2))))/2
    distSatellite2User = UC*cosd(angleOfElevation) - sqrt(((UC*cosd(angleOfElevation)).^2) - ...
                            ((UC.^2) - (SC.^2)));

    % Angle between SC and SU
    angleBetnSCnSU = acosd(((SC.^2) + (distSatellite2User.^2) - (UC.^2))/(2*(SC.*distSatellite2User)));
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
    satelliteAngletoNorth = acosd(sum(velocityVector.*northAxisVector)/(sqrt((sum(velocityVector.^2))*(sum(northAxisVector.^2)))));
end

%% Function Details
% Name              - computeDoppler
% Details           - Computes Doppler
% Input Parameters  - sysConfig         - Structure with all the parameters
%                                         needed
% Output Parameters - observedDoppler   - Observed Doppler
function observedDoppler = computeDoppler(sysConfig)

    % If using Genie based Doppler Computation
    if sysConfig.genieBasedDopplerComp
        % Identifying the plane of Satellite (S), Centre of Earth (C), and User (U)
        % Converting Latitudes and Longitudes to Cartesian Coordinates
        earthCentreCoordinate = [0, 0, 0];
        satelliteCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
        ueCoordinates = convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation);

        % Given coordinates, plane being Ax+By+Cz+D=0
        % vector v1 = satelliteCoordinates - earthCentreCoordinate, v2 = ueCoordinates - earthCentreCoordinate
        % (A, B, C) = v1 x v2 (x => Cross Product)
        v1 = satelliteCoordinates - earthCentreCoordinate;
        v2 = ueCoordinates - earthCentreCoordinate;
        A = v1(2)*v2(3) - v2(2)*v1(3);
        B = v1(3)*v2(1) - v2(3)*v1(1);
        C = v1(1)*v2(2) - v2(1)*v1(2);
        % Since the plane should pass through Origin (defined from earthCentreCoordinate), D = 0
        % D = -A*earthCentreCoordinate(1) - B*earthCentreCoordinate(2) - C*earthCentreCoordinate(3);

        % Projection angle of velocity vector and CUS Plane
        velocityVector = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc2) - ...
                            convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
        cusPlaneVector = [A, B, C];
        sysConfig.velocityProjectionAngleOntoCusPlaneGenie = asind(sum(velocityVector.*cusPlaneVector)/(sqrt(...
                                                    (sum(velocityVector.^2))*(sum(cusPlaneVector.^2)))));

        % Projection of Velocity Vector onto CUS Plane
        velocityCusPlaneProjection = velocityVector - ...
                            cusPlaneVector.*(sum(velocityVector.*cusPlaneVector)/sum(cusPlaneVector.^2));

        % Angle between Satellite, User and Projected Velocity Vector
        pathVector = ueCoordinates - satelliteCoordinates;
        sysConfig.angleOfDepressionGenie = acosd(sum(pathVector.*velocityCusPlaneProjection)/(sqrt(...
                                    (sum(pathVector.^2))*(sum(velocityCusPlaneProjection.^2)))));
        velocityOnPathProjectionVector = pathVector.*(sum(velocityCusPlaneProjection.*pathVector)/sum(pathVector.^2));

        % Observed Doppler
        observedDoppler = sqrt(sum(velocityOnPathProjectionVector.^2));

    else

        % Practical Doppler Computation
        observedDoppler = 20; % Dummy Code

    end
end

%% Function Details
% Name              - computeAoA
% Details           - Computes Angle of Arrival
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - observedAoA    - Angle of Arrival
function observedAoA = computeAoA(sysConfig)

    % Genie AoA
    if sysConfig.genieBasedAoAComp
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
    else

        % Practical AoA Estimation
        observedAoA = 30; % Dummy Code

    end
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