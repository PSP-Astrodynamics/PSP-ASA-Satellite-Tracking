function stateDerivative = satelliteDynamics(JulianDate, satelliteState, Spacecraft, Environment)
%=========================================================================%
% FUNCTION: satelliteDynamics
% AUTHOR: Your Name (your email), 2021
% DESCRIPTION: Given a satellite's position and velocity at a certain time, 
% this function computes the velocity and acceleration of the satellite.
% INPUTS:
%   > JulianDate = Julian date at the current time (including the fraction
%     of the Julian date)
%   > satelliteState = a 6x1 vector containing the position vector of the
%     satellite in the J2000 Earth Mean Equator reference frame [km] for
%     the first 3 components, and the velocity in the same frame [km/s] for
%     the last 3 components.
%   > Spacecraft: a spacecraft object containing the following properties
%               .mass = spacecraft mass [kg]
%               .Aref = spacecraft reference area [m^2]
%               .Cd   = spacecraft drag coefficient [unitless]
%               .JDepoch = Julian date of initial epoch
%               .position = spacecraft cartesian position at epoch in the
%               J2000 EME frame [km]
%               .velocity = spacecraft cartesian velocity at epoch in the
%               J2000 EME frame [km/sec]
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
% OUTPUTS:
%   > stateDerivative = a 6x1 vector containing the velocity in the J2000 
%     Earth Mean Equator reference frame [km/s] for the first 3 components
%     and the acceleration [km/s^2] in the same frame as the last 3 
%     components.
%=========================================================================%

%% Unpack variables
r = satelliteState(1:3); % Position [km]
v = satelliteState(4:6); % Velocity [km/s]

%% Compute two-body gravity acceleration
accelTwoBody = -Environment.EarthMU/norm(r)^3*r;

%% Compute J2 gravity perturbation
J2coeff = -3/2*Environment.EarthJ2*Environment.EarthMU*Environment.EarthEquatorialRadius^2/norm(r)^5;
accelJ2(1,1) = J2coeff*(1-5*(r(3)/norm(r))^2)*r(1);
accelJ2(2,1) = J2coeff*(1-5*(r(3)/norm(r))^2)*r(2);
accelJ2(3,1) = J2coeff*(3-5*(r(3)/norm(r))^2)*r(3);

%% Compute aerodynamic drag force

% Compute Earth-fixed position and velocity of the satellite   
[recef,vecef,DCM]  = convertECI2ECEF(r,v,JulianDate,Environment);
vAtmos = DCM'*vecef; % Atmospheric relative velocity in ECI frame

% Compute latitude, longitude, and altitude of the satellite
% NOTE: longitude will vary between -180 and 180 degrees in accordance to
% what is expected by the NRL MSIS model. 
latitude=real(asind(recef(3)/norm(recef))); %[deg]
longitude = atan2(recef(2),recef(1))*180/pi; %[deg]
altitude = (norm(r)-Environment.EarthEquatorialRadius)*1000; %[meters]!!!
if altitude < 100e3
    error('Satellite entered atmosphere')
end

% Convert Julian date to Gregorian date and extract, year, dayOfYear, and
% UTC seconds
GregorianDateTime = datetime(2000,1,1,0,0,0) + days(JulianDate - 2451544.5);
yr = year(GregorianDateTime);
dayOfYear = day(GregorianDateTime,'dayofyear');
UTCseconds = (JulianDate + 0.5 - floor(JulianDate + 0.5))*86400;

% Compute atmospheric density
%if altitude > 1e6
%    airDen = 0;
%else
[~, rho] = atmosnrlmsise00(altitude,latitude,longitude,yr,dayOfYear,UTCseconds, Environment.f107Average, Environment.f107Daily, Environment.magneticIndex);
airDen = rho(6); % Total mass density [kg/m^3]
%end

% Compute atmospheric drag acceleration
accelDrag = 0.5*airDen*norm(vAtmos)^2 * Spacecraft.Cd*Spacecraft.Aref / Spacecraft.mass * (-vAtmos/norm(vAtmos))*1000;

%% Compile output
accelTotal = accelTwoBody + accelJ2 + accelDrag;
stateDerivative = [v;accelTotal];

end
