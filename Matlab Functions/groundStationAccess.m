function [satElevation,satRange,satDoppler] = groundStationAccess(groundStation,propJulianDates,propStates,Environment)
%=========================================================================%
% FUNCTION: groundStationAcess
% AUTHOR: Your Name (your email), 2021
% DESCRIPTION: given a ground station and orbital ephemerides, computes
% the satellite's elevation above the horizon of the ground station, it's
% range, and the Doppler shift.
% INPUTS:
%   > groundStation = ground station object containing properties:
%                  .latitude = latitude of ground station [deg]
%                  .longitude = longitude of ground station [deg]
%                  .altitude = altitude of ground station [km]
%   > propJulianDates = 1xN Julian date vector corresponding to the orbit
%   ephemerides.
%   > propStates = 6xN vector containing the position vector of the
%     satellite in the J2000 Earth Mean Equator reference frame [km] for
%     the first 3 components, and the velocity in the same frame [km/s] for
%     the last 3 components.
%   > Environment = environment object containing the following properties
%               .EarthMU = Earth standard gravitational parameter
%               [km^3/s^2]
%               .EarthJ2 = Earth J2 gravity constant [unitless]
%               .EarthPolarRadius = Earth polar radius [km]
%               .EarthEquatorialRadius = Earth equatorial radius [km]
%               .EarthRotationRate = Earth rotational rate [rad/s]
%               .f107Daily = Solar F10.7 cm radio flux [SFU]
%               .f107Average = 81-day average F10.7 cm radio flux [SFU]
%               .magneticIndex = Geomagnetic activity index [unitless]
% OUPUTS:
%   > satElevation = 1xN vector of the elevation angle of the satellite.
%   > satRange = 1xN vector of the range to the satellite [km]
%   > satDoppler = 1xN vector of the Doppler shift fraction (frequency
%   observed divided by source frequency)
%=========================================================================%

%% Unpack variables
r = propStates(1:3,:); % Position in ECI frame [km]
v = propStates(4:6,:); % Velocity in ECI frame [km/s]
lat = groundStation.latitude; % Latitude [deg]
lon = groundStation.longitude; % Longitude [deg]
alt = groundStation.altitude; % Altitude [km]
Requator = Environment.EarthEquatorialRadius; % [km]
Rpolar = Environment.EarthPolarRadius; % [km]

%% Ground station position vector and normal vector in Earth-fixed frame

% Earth Ellipsoid Flattening 
f = (Requator - Rpolar)/Requator; % Eq. 5.77 in reference

% Radius of curvature of in the prime vertical
N = Environment.EarthEquatorialRadius/sqrt(1 - f*(2 - f)*sind(lat)^2); % Eq. 5.84 in reference

% Cartesian location of ground station (Eq. 5.83 in reference)
rGS = [(N+alt)*cosd(lat)*cosd(lon); 
       (N+alt)*cosd(lat)*sind(lon); 
       ((1-f)^2*N+alt)*sind(lat) ];

% Normal vector
%nHorizon = 2*[rGS(1)/Requator^2; rGS(2)/Requator^2; rGS(3)/Rpolar^2];
nHorizon = rGS/norm(rGS);
   
%% Satellite positions and velocities in the Earth-fixed frame
for ii = 1:numel(propJulianDates)
    [recef(:,ii),vecef(:,ii),~]  = convertECI2ECEF(r(:,ii),v(:,ii),propJulianDates(ii),Environment);
end

%% Determine elevation, range, and Doppler shift
for ii = 1:numel(propJulianDates)
    
    % Elevation
    rStation2Sat = recef(:,ii) - rGS; % Vector from station to satellite
    satElevation(ii) = 90 - acosd(dot(rStation2Sat,nHorizon)/(norm(rStation2Sat)*norm(nHorizon)));
    
    % Range
    satRange(ii) = norm(rStation2Sat);
    
    % Doppler
    c = 2.998e5; % Speed of light [km/s]
    relativeVelocity = dot(vecef(:,ii),rStation2Sat/norm(rStation2Sat));
    satDoppler(ii) = c/(c + relativeVelocity);
end


end