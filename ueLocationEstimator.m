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
sysConfig.satelliteLoc2           = [38.50189, -121.520728]; % [Latitude, Longitude]

% Variables regarding UE
% Radius of Earth
sysConfig.earthRadius             = 6400; % m
% User Location
sysConfig.ueLocation              = [39.709380, -121.100728];% [39.195778, -123.280982];% [39.034394, -119.280038]; % [37.389108, -122.143192]; % [Latitude, Longitude] 
% User Azimuth wrt North
sysConfig.ueAzimuthWN             = 30; % Degrees towards East

% Carrier Frequency
sysConfig.carrierFreq             = 11.72e9; % Hz

% Synchronization
sysConfig.synchronized            = false;
sysConfig.errorBiasInToA          = 0.0; % Percentage

% Movement Sign of SNR (shall be reassigned to meet the requirements)
% +1 -> Satellite Going Away or SNR historical trend for this satellite is dereasing.
% -1 -> Satellite Coming Towards or SNR historical trend for this satellite is increasing.
sysConfig.movementSign            = +1;

% AoA Computation
sysConfig.genieBasedAoAComp       = true;
sysConfig.errorBiasInAoA          = 0.0; % Degree

% Doppler Computation
sysConfig.genieBasedDopplerComp   = true;
sysConfig.errorPercInDoppler      = 0.0; % Percentage
sysConfig.maxDopplerError         = 600; % Hz

% Loop on Dopplers
snrArray = 0 :2/3 : 20;
% consDopplerPercChanges = 0.0:0.0002:0.006;
% consAoABiasChanges = 0.0:0.0033:0.2;
consDopplerPercChanges = [sqrt(2e-3):(sqrt(6e-4)-sqrt(2e-3))/7:sqrt(6e-4), ...
                            sqrt(6e-4)+((sqrt(5.5e-4)-sqrt(6e-4))/8):(sqrt(5.5e-4)-sqrt(6e-4))/8:sqrt(5.5e-4), ...
                            sqrt(5.5e-4)+((sqrt(5e-4)-sqrt(5.5e-4))/15):(sqrt(5e-4)-sqrt(5.5e-4))/15:sqrt(5e-4)]./20;
consAoABiasChanges = [0.18:(0.157-0.18)/7:0.157, ...
                        0.157+((0.125-0.157)/8):(0.125-0.157)/8:0.125, ...
                        0.125+((0.119-0.125)/8):(0.119-0.125)/8:0.119, ...
                        0.119+((0.116-0.119)/7):(0.116-0.119)/7:0.116]/1.3;
mseDopplers = zeros(1, length(consDopplerPercChanges));
meanDopplers = zeros(1, length(consDopplerPercChanges));
for currDopplerPercChangeIdx = 1 : length(consDopplerPercChanges)
% for currDopplerPercChangeIdx = 1 : length(consAoABiasChanges)

    sysConfig.errorPercInDoppler = consDopplerPercChanges(currDopplerPercChangeIdx);
    sysConfig.errorBiasInAoA = consAoABiasChanges(currDopplerPercChangeIdx);

    %% Looping on Simulations
    numSimulations = 100000;
    mseErrorSet = zeros(1, numSimulations);
    dopplerSet = zeros(1, numSimulations);
    mErrDoppler = zeros(1, numSimulations);
    for simIdx = 1 : numSimulations
        % UE Location
        sysConfig.ueLocation              = [32 + 10*rand(1), -124 + 4*rand(1)];
    
        % Calculating AoA
        sysConfig = computeAoA(sysConfig);
        
        % Calculating Doppler
        sysConfig = computeDoppler(sysConfig);
        
        % Fetching UE Location
        sysConfig = getUeLoc(sysConfig);
    
        % MSE Pooling
        mseErrorSet(simIdx) = sqrt(sum((convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation) - ...
                    convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.estimatedUeLocation)).^2));
        
        dopplerSet(simIdx) = sysConfig.observedDopplerGenie;
        mErrDoppler(simIdx) = sysConfig.errorDoppler;
    
        % % Dispalying Results and Accuracy 
        % [sysConfig.ueLocation, sysConfig.estimatedUeLocation, ...
        %     sqrt(sum((convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation) - ...
        %                 convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.estimatedUeLocation)).^2))]

    end % End of loop on Simulations
    mseDopplers(currDopplerPercChangeIdx) = mean(mseErrorSet);
    meanDopplers(currDopplerPercChangeIdx) = mean(mErrDoppler);

end

% plot(consDopplerPercChanges, mseDopplers, "LineWidth", 2);
% data = sortrows([meanDopplers.', mseDopplers.']);
% plot(data(:,1), data(:,2), "LineWidth", 2);
% grid on;
% xlabel("Doppler Error (Hz)");
% title("Position Error vs. Doppler Error");
% xlabel("Doppler % Change for errors");
% title("Position Error vs. Doppler % Error");
% plot(consAoABiasChanges, mseDopplers, "LineWidth", 2);
% grid on;
% xlabel("AoA Change for errors");
% title("Position Error vs. AoA Error");
plot(snrArray, mseDopplers, "LineWidth", 2);
grid on;
xlabel("SNR (dB)");
title("Position Error vs. SNR (dB)");
ylabel("Mean Position Error (Km)");

%% Function Details
% Name              - getUeLoc
% Details        s   - Computes User Location using Iterative procedure
% Input Parameters  - sysConfig         - Structure with all the parameters
%                                         needed
% Output Parameters - ueLocation        - User Location
function ueLocation = getUeLoc(sysConfig)

    % Achieving the distance of the Satellite(S) to User(U)
    [sysConfig.distSatellite2User, sysConfig.angleBetnSCnSU] = getDistUser2Satellite(sysConfig);

    % Angle between velocity vector projected on SCU Plane 
    sysConfig.phi = 90 + sysConfig.movementSign*sysConfig.angleBetnSCnSU;

    % Projection Angle of Satellite Velocity and Plane containing User,
    % Satellite and Centre of Earth
    sysConfig.satelliteProjectionAngle = acosd(max(min(sysConfig.observedDoppler*3e8/(sysConfig.satelliteVelocity*sysConfig.carrierFreq*cosd(sysConfig.phi)), 1), -1));

    % Get the UE Location
    sysConfig.newUeLoc = computeUeCoordinates(sysConfig);

    % UE Location
    ueLocation = sysConfig.newUeLoc;
end


%% Function Details
% Name              - computeUeCoordinates
% Details           - Computes UE Location using it's distance, angle of
%                     depression wrt satellite and user's angle of arrival
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                       needed
% Output Parameters - ueLocation     - Structure with all the parameters
%                                         needed with User's Location Coordinates
function sysConfig = computeUeCoordinates(sysConfig)
    
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
    velocityVector = velocityVector./sqrt(sum(velocityVector.^2));
    tangentialVelocityVecPlane = convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorStart);
    velocityPerpVectorInCusPlane = cross(tangentialVelocityVecPlane, velocityVector);
    if sysConfig.satelliteProjectionAngle > 90
        sysConfig.satelliteProjectionAngle = sysConfig.satelliteProjectionAngle - 180;
    end

    velocityVectorProjectionInCusPlane = velocityVector*cosd(leftRightSign*sysConfig.satelliteProjectionAngle) + ...
                    sind(leftRightSign*sysConfig.satelliteProjectionAngle)*...
                    velocityPerpVectorInCusPlane./sqrt(sum(velocityPerpVectorInCusPlane.^2));

    % Direction on projected velocity in the path direction
    satelliteCentreVector = tangentialVelocityVecPlane;
    orthoCusPlaneVector = cross(satelliteCentreVector, velocityVectorProjectionInCusPlane);
    orthoVelocityInCusPlaneVector = cross(orthoCusPlaneVector, velocityVectorProjectionInCusPlane);
    pathDirection = velocityVectorProjectionInCusPlane*cosd(sysConfig.phi) + ...
                        sind(sysConfig.phi)*orthoVelocityInCusPlaneVector./sqrt(sum(orthoVelocityInCusPlaneVector.^2));
    pathDirectionUnitNorm = pathDirection./sqrt(sum(pathDirection.^2));

    % Finding the UE Location in the path direction from Satellite Location
    ueLocationCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, velocityVectorStart) + ...
                               pathDirectionUnitNorm*sysConfig.distSatellite2User;

    % Converting UE Location to Latitude and Longitude Co-ordinates
    
    % [convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation), ueLocationCoordinates]
    % [sqrt(sum((convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation) - ...
    %             ueLocationCoordinates).^2))]

    sysConfig.estimatedUeLocation = convCartesian2LatLong(sysConfig.earthRadius, ueLocationCoordinates);
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

    if sysConfig.synchronized
        % Distance Between User and Satellite
        % is directly related to timeOfArrival(usec) due to LOS
        distSatellite2User = sysConfig.timeOfArrival*3e8/1e6/1e3; % km
    else
        % Angle of Elevation
        angleOfElevation = abs(sysConfig.observedAoA) + 90;

        % We know SC, UC, and angle between UC, SU (angleOfElevation)
        % SU can be solved using the equation solving:
        % SU^2 - 2*SU*UC*cos(angleOfElevation) + UC^2 - SC^2 = 0
        % SU = (2*UC*cosd(angleOfElevation) - sqrt(((2*UC*cosd(angleOfElevation))^2) - ...
        %           4*((UC^2)-(SC^2))))/2
        distSatellite2User = abs(UC*cosd(angleOfElevation) + sqrt(((UC*cosd(angleOfElevation)).^2) - ...
                                ((UC.^2) - (SC.^2))));
    end

    % Angle between SC and SU
    angleBetnSCnSURaw = ((SC.^2) + (distSatellite2User.^2) - (UC.^2))/(2*(SC.*distSatellite2User));
    angleBetnSCnSU = acosd(max(min(angleBetnSCnSURaw, 1), -1));

    % Update the distance for error resilience
    if abs(angleBetnSCnSURaw - cosd(angleBetnSCnSU)) > 1e-3
        distSatellite2User = SC*cosd(angleBetnSCnSU) + sqrt((UC.^2) - ((SC*cosd(angleBetnSCnSU)).^2));
    end

end

%% Function Details
% Name              - computeDoppler
% Details           - Computes Doppler
% Input Parameters  - sysConfig         - Structure with all the parameters
%                                         needed
% Output Parameters - sysConfig         - Structure with all the parameters
%                                         needed with Observed Doppler, 
%                                         additional genie parameters for
%                                         verification
function sysConfig = computeDoppler(sysConfig)

    % If using Genie based Doppler Computation
    if sysConfig.genieBasedDopplerComp
        % Identifying the plane of Satellite (S), Centre of Earth (C), and User (U)
        % Converting Latitudes and Longitudes to Cartesian Coordinates
        earthCentreCoordinate = [0, 0, 0];
        satelliteCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
        ueCoordinates = convLatLong2Cartesian(sysConfig.earthRadius, sysConfig.ueLocation);

        % Given coordinates, plane being Ax+By+Cz+D=0
        % vector v1 = satelliteCoordinates - earthCentreCoordinate, v2 = ueCoordinates - earthCentreCoordinate
        v1 = satelliteCoordinates - earthCentreCoordinate;
        v2 = ueCoordinates - earthCentreCoordinate;
        cusPlaneNormVector = cross(v1, v2);
        % Since the plane should pass through Origin (defined from earthCentreCoordinate), D = 0
        % D = -A*earthCentreCoordinate(1) - B*earthCentreCoordinate(2) - C*earthCentreCoordinate(3);

        % Projection angle of velocity vector and CUS Plane
        velocityVector = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc2) - ...
                            convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc1);
        velocityVector = sysConfig.satelliteVelocity*velocityVector/sqrt(sum(velocityVector.^2));
        sysConfig.velocityProjectionAngleOntoCusPlaneGenie = asind(sum(velocityVector.*cusPlaneNormVector)/(sqrt(...
                                                    (sum(velocityVector.^2))*(sum(cusPlaneNormVector.^2)))));

        if sysConfig.velocityProjectionAngleOntoCusPlaneGenie < 0
            sysConfig.velocityProjectionAngleOntoCusPlaneGenie = ...
                sysConfig.velocityProjectionAngleOntoCusPlaneGenie + 180;
        end

        % Projection of Velocity Vector onto CUS Plane
        velocityCusPlaneProjection = velocityVector - ...
                            cusPlaneNormVector.*(sum(velocityVector.*cusPlaneNormVector)/sum(cusPlaneNormVector.^2));

        % Angle between Satellite, User and Projected Velocity Vector
        pathVector = ueCoordinates - satelliteCoordinates;
        sysConfig.angleOfDepressionGenie = acosd(sum(pathVector.*velocityCusPlaneProjection)/(sqrt(...
                                    (sum(pathVector.^2))*(sum(velocityCusPlaneProjection.^2)))));
        if sysConfig.angleOfDepressionGenie > 90
            sysConfig.movementSign = +1;
        else
            sysConfig.movementSign = -1;
        end
        velocityOnPathProjectionVector = pathVector.*(sum(velocityCusPlaneProjection.*pathVector)/sum(pathVector.^2));

        % Observed Doppler
        sysConfig.observedDopplerGenie = sqrt(sum(velocityOnPathProjectionVector.^2))*sysConfig.carrierFreq/3e8;

        % Error Injection
        sysConfig.errorDoppler = min(sysConfig.errorPercInDoppler*sysConfig.observedDopplerGenie, sysConfig.maxDopplerError);
        sysConfig.observedDoppler = (1-2*randi([0, 1]))*sysConfig.errorDoppler + ...
                                        sysConfig.observedDopplerGenie;
    else

        % Practical Doppler Computation
        sysConfig.observedDoppler = 20; % Dummy Code

    end
end

%% Function Details
% Name              - computeAoA
% Details           - Computes Angle of Arrival
% Input Parameters  - sysConfig      - Structure with all the parameters
%                                      needed
% Output Parameters - sysConfig      - Structure with all the parameters
%                                      needed with Observed AoA, 
%                                      movement sign
function sysConfig = computeAoA(sysConfig)

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

        % Side of UE
        nextSatelliteCoordinates = convLatLong2Cartesian(sysConfig.distSatellite2Centre, sysConfig.satelliteLoc2);
        velocityVector = nextSatelliteCoordinates - satelliteCoordinates;
        unitVelocityVector = velocityVector./sqrt(sum(velocityVector.^2));

        satelliteEarthCentreVector = -satelliteCoordinates./sqrt(sum(satelliteCoordinates.^2));
        planeOfSatelliteVelocityAndCentre = cross(unitVelocityVector, satelliteEarthCentreVector);
        sideOfUe = sign(sum(planeOfSatelliteVelocityAndCentre.*ueCoordinates));

        % Time of Arrival
        sysConfig.timeOfArrivalGenie = SU*1e9/3e8; % usec

        % Angle of Arrival
        sysConfig.observedAoAGenie = sideOfUe*(acosd(((UC.^2) + (SU.^2) - (SC.^2))./(2*(UC.*SU))) - 90);

        % Error Injection
        sysConfig.observedAoA = (1-2*randi([0, 1]))*sysConfig.errorBiasInAoA + ...
                                        sysConfig.observedAoAGenie;
        sysConfig.timeOfArrival = (1-2*randi([0, 1]))*sysConfig.errorBiasInToA + ...
                                        sysConfig.timeOfArrivalGenie;
    else

        % Practical AoA Estimation
        sysConfig.observedAoA = 30; % Dummy Code

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
    latitude = asind(max(min(cartesianCoordinates(3)/radius, 1), -1));
    currSign = sign(asind(cartesianCoordinates(2)/radius/cosd(latitude)));
    longitude = currSign*acosd(max(min(cartesianCoordinates(1)/radius/cosd(latitude), 1), -1));

    % Putting them together
    latLongPair = [latitude, longitude];
end