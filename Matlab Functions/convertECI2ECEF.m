function [recef,vecef,DCM]  = convertECI2ECEF(r,v,JulianDate,Environment)
%=========================================================================%
% FUNCTION: convertECI2ECEF
% AUTHOR: Your Name (your email), 2021
% DESCRIPTION: Computes the direction cosine matrix (DCM) to rotate a
% column vector into the Earth-fixed frame and applies this to a position
% and velocity vector.
% INPUTS:
%   > r = 3x1 position vector in the J2000 EME frame
%   > v = 3x1 velocity vector in the J2000 EME frame
%   > JulianDate = Julian date at the current time (including the fraction
%     of the Julian date)
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
%   > recef = 3x1 position vector in the Earth fixed frame
%   > vecef = 3x1 velocity vector in the Earth fixed frame
%   > DCM   = direction cosine matrix to rotate a  a
%     column vector from the J2000 frame into the Earth-fixed frame
%=========================================================================%

% Compute Julian centuries from 2000 (T) and UT1, assuming UTC = UT1
T = (JulianDate-2451545.0)/36525;
UT1 = (JulianDate + 0.5 - floor(JulianDate + 0.5))*86400;

% Compute the Greenwich Mean Sidereal Time (GMST)
GreenwichMeanSiderealTime = 24110.54841 + 8640184.812866*T + ...
                    0.093104*T^2 - 0.0000062*T^3 + 1.002737909350795*UT1;
GreenwichMeanSiderealTime = mod(GreenwichMeanSiderealTime,86400);

% Express GMST as an hour angle between 0 and 2pi radians
GreenwichMeanSiderealTime = GreenwichMeanSiderealTime*2*pi/86400;  %(rad)

% Direction cosine matrix - simple rotation about z axis
DCM = [ cos(GreenwichMeanSiderealTime), sin(GreenwichMeanSiderealTime),0;...
       -sin(GreenwichMeanSiderealTime), cos(GreenwichMeanSiderealTime),0;...
                   0,                                 0,               1];

% Compute Earth-fixed position and velocity of the satellite               
recef = DCM*r;
vecef = DCM*(v + cross(r,Environment.EarthRotationRate*[0;0;1]));

end