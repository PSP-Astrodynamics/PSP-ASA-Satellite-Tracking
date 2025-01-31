function [propJulianDates, propStates] = propagateSatellite(tspan, Spacecraft, Environment)
%=========================================================================%
% FUNCTION: propagateSatellite.m
% AUTHOR: Your Name (your email), 2021
% DESCRIPTION: Integrates the position of a spacecraft using a simple
% dynamics model.
% INPUT:
%   > tspan: a vector of times [in sec] after the initial spacecraft epoch
%   at which the propagated satellite state is desired.
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
% OUTPUT:
%   > propJulianDates: 1xN vector of Julian dates for the output times
%   > propStates: 6xN vector of propagated states in the J2000 EME frame.
%     The first three components are the position [km] and the next three
%     components are the velocity [km/sec]
%=========================================================================%

%% Propagation Settings
maxStep = 60; % Maximum step size [sec]

%% Configure Propagator
odeOptions = odeset('MaxStep',maxStep);
state0 = [Spacecraft.position; Spacecraft.velocity];
odefun = @(t,satelliteState) satelliteODE(t, satelliteState);

%% Run Propagator
[propSecs,propStates] = ode113(odefun,tspan,state0,odeOptions);
propJulianDates = Spacecraft.JDepoch + propSecs/86400; 

%% Propagation ODE function
function stateDerivative = satelliteODE(t, satelliteState)

JulianDate = Spacecraft.JDepoch + t/86400;
stateDerivative = satelliteDynamics(JulianDate, satelliteState, Spacecraft, Environment);

end

end
