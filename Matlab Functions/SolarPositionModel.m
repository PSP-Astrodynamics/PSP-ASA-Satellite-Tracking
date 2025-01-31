function rSun = SolarPositionModel(JD)
%=========================================================================%
% FUNCTION: SolarPositionModel.m
% AUTHOR: Justin Mansell (jmansell@purdue.edu), 2021
% DESCRIPTION: Given a Julian date, returns the position of the Sun in the
% J2000 Earth Mean Equator inertial frame in units of km.
%=========================================================================%

% Compute Julian centuries from 2000 (T) and UT1, assuming UTC = UT1
T = (JD-2451545.0)/36525;

% Determine angle of Sun along ecliptic
OmegaPlusOmega = 282.9400; % Ω + w [deg]
meanAnom = 357.5256+35999.049*T; % Mean anomaly [deg]
lambdaSun = OmegaPlusOmega + meanAnom + 6892/3600*sind(meanAnom) + ...
                                        72/3600 + sind(2*meanAnom); %[deg]

% Determine distance to Sun
distSun = (149.619 - 2.499*cosd(meanAnom) - 0.021*cosd(2*meanAnom))*1e6; %[km]

% Determine position vector
obliquityOfEcliptic = 23.43929111; % Obliquity of the ecliptic [deg]
rSun = [distSun.*cosd(lambdaSun);
        distSun.*sind(lambdaSun)*cosd(obliquityOfEcliptic);
        distSun.*sind(lambdaSun)*sind(obliquityOfEcliptic)];

end