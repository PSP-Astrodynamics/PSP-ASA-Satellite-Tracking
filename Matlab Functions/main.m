%=========================================================================%
% PROGRAM: main.m
% AUTHOR: Instructors, 2023
% DESCRIPTION: This file provides orbit propagation and related info 
% for the Satellite Tracking Project. 
%=========================================================================%
clear; close all
%--------------------------USER SETTINGS----------------------------------%

%% Define Satellite TLE - UPDATE THESE VALUES
tleLine1 = '1 51085U 22002DF  25032.54719027  .00048242  00000+0  10702-2 0  9993';
tleLine2 = '2 51085  97.3538 102.9956 0006702 315.2737  44.7970 15.43560689169443';

%% Define Satellite Properties - UPDATE THESE VALUES
Spacecraft.mass = 2;      % Satellite mass [kg]
Spacecraft.Aref = .03405; % Satellite reference area [m^2]
Spacecraft.Cd   = 0.2;    % Satellite drag coefficient [unitless]

%% Define Ground Station Properties - UPDATE THESE VALUES
Station.latitude     = 50.359;    % Ground station latitude [deg]
Station.longitude    = 30.388;    % Ground station longitude [deg]
Station.altitude     = 0.2;       % Ground station altitude [km]
Station.minElevation = 5;         % Minimum elevation for acquisition [deg]
Station.freq         = 435000000; % Station frequency [Hz]

%------------------------END USER SETTINGS--------------------------------%

%% SATELLITE POSITION/VELOCITY CALCULATIONS
[rxyz, velxyz, alt, JD] = getInitialStateVectorFunc(tleLine1, tleLine2);
Spacecraft.JDepoch      = JD;     % Satellite Julian date at epoch
Spacecraft.position     = rxyz;   % Satellite initial ECI position [km]
Spacecraft.velocity     = velxyz; % Satellite initial ECI velocity [km/s]

%% ENVIRONMENT PROPERTIES
[Environment.f107Daily, Environment.f107Average, Environment.magneticIndex] = getAverageF107(); 
% [Solar F10.7 cm radio flux [SFU], 81-day average F10.7 cm radio flux [SFU], Geomagnetic activity index [unitless]]
Environment.EarthMU = 3.9860044189e5; % Earth gravitational parameter [km^3/s^2]
Environment.EarthJ2 = 1082.63e-6; % J2 oblateness parameter [unitless]
Environment.EarthPolarRadius = 6356.752; % [km]
Environment.EarthEquatorialRadius = 6378.1363; % [km]
Environment.EarthRotationRate = 7.2921159e-5; % Earth's rotation rate [rad/s]

%% Propagate Properties
outputTimeStep = 10; % Time step at which to output data [sec]
maxTime = 86400*3;  % Maximum time to propagate to [sec], e.g., 3 days

%% Run Propagation
tspan = 0:outputTimeStep:maxTime; % Propagation time vector
[propJulianDates, propStates] = propagateSatellite(tspan, Spacecraft, Environment);
propStates = propStates'; % Propagated state vector

%% Compute ground contacts
[satElevation,satRange,satDoppler] = groundStationAccess(Station,propJulianDates,propStates,Environment);
% Find at which times the satellite is visible to ground station
visibleInds = satElevation > Station.minElevation;
visibleJulianDates = propJulianDates(visibleInds);

% Compute range and Doppler frequency where satellite is visible to station
rangeObs = satRange(visibleInds);
dopplerObs = satDoppler(visibleInds)*Station.freq;

%% Determine eclipse
for ii=1:numel(propJulianDates)
    JD = propJulianDates(ii); % Current Julian date
    rSat = propStates(1:3,ii); % Current satellite position vector
    rSun = SolarPositionModel(JD); % Solar position vector
    inEclipse(ii) = eclipseCheck(rSun,rSat);
end
inEclipse = inEclipse';

%% Plots

% Plot elevation above ground station
t = days(propJulianDates - 2451544.5) + datetime(2000,1,1,0,0,0);
figure(1)
plot(t,satElevation,'k-')
ylabel('Elevation [deg]')
title('Elevation Above Ground Station Horizon')
xlim([t(1) t(end)])
hold on
yline(Station.minElevation);
grid on

%% Plot passes with contrained elevation
passtimeUp = t(1);
passtimeDown = t(1);
passtimePeakTimes = [];
passtimePeakValues = [];

for i = 1:(length(satElevation)-1)
    if (satElevation(i) <= Station.minElevation) && (satElevation(i+1) >= Station.minElevation)
        passtimeUp = [passtimeUp; t(i)];
    end
    if (satElevation(i) >= Station.minElevation) && (satElevation(i+1) <= Station.minElevation)
        passtimeDown = [passtimeDown; t(i)];
    end
end

passtime = [passtimeUp(2:end), passtimeDown(2:end)];

% Record the peak elevation between passtimeUp and passtimeDown for each pass
for j = 1:length(passtime(:,1))
    % Find indices for the current pass
    idxStart = find(t >= passtime(j,1), 1, 'first');
    idxEnd = find(t <= passtime(j,2), 1, 'last');
    
    % Extract elevation values for this pass
    [peakValue, peakIdx] = max(satElevation(idxStart:idxEnd));
    
    % Record peak elevation value and corresponding time
    passtimePeakValues = [passtimePeakValues; peakValue];
    passtimePeakTimes = [passtimePeakTimes; t(idxStart + peakIdx - 1)];
end
 
 

% Plot range and Doppler shift
t2 = days(visibleJulianDates - 2451544.5) + datetime(2000,1,1,0,0,0);
figure(2)
subplot(2,1,1)
plot(t2,rangeObs,'ko')
ylabel('Range [km]')
grid on
title('Range and Recieved Frequency')
subplot(2,1,2)
plot(t2,dopplerObs/1e6,'ko')
ylabel('Frequency [MHz]')
grid on

% Plot eclipse state
figure(3)
plot(t,inEclipse,'k')
yticks([0 1])
yticklabels({'Sunlit','Eclipse'})
xlim([t(1) t(end)])

% Plot geocentric altitude
figure(4)
hG = vecnorm(propStates(1:3,:),2,1) - Environment.EarthEquatorialRadius;
plot(t,hG,'k')
xlim([t(1) t(end)])
ylabel('Geocentric Altitude [km]')
grid on

% Determine state at query time
tq = datetime(2024,11,28,15,5,0);
%tq = datetime(2460639.333973, 'ConvertFrom', 'juliandate');
positionEstimate = [interp1(t,propStates(1,:),tq);
                    interp1(t,propStates(2,:),tq);
                    interp1(t,propStates(3,:),tq)];
velocityEstimate = [interp1(t,propStates(4,:),tq);
                    interp1(t,propStates(5,:),tq);
                    interp1(t,propStates(6,:),tq)];
h_query = norm(positionEstimate)-6378.1363;

%% 3d propagation
figure
plotStaticPropagation(propStates')
hold on


%% RAAN
R_E = 6378.1363; % km
mu_E =  398600.4415; % km^3/s^2
[a, ecc, incl, RAAN, omega, theta_star, theta] = deal(zeros(1,ii));
for i=1:ii
    % [a(i), ecc(i), incl(i), Omega(i), omega(i), theta_star(i), trueLon(i), theta(i), lonPer(i)] =...
    %     ijk2keplerian(r_ijk, v_ijk);
    [a(i), ecc(i), incl(i), RAAN(i), omega(i), theta_star(i), ~, theta(i), ~] =...
        ijk2keplerian(1000*propStates(1:3,i), 1000*propStates(4:6,i));
    a(i) = a(i)/1000; % km

    r(1:3,i) = propStates(1:3,i); % Position [km]
    v(1:3,i) = propStates(4:6,i); % Velocity [km/s]
    J2coeff(i) = -3/2*Environment.EarthJ2*Environment.EarthMU*Environment.EarthEquatorialRadius^2/norm(r)^5;
    
    IP(i) = 2*pi*sqrt((a(i)^3)/mu_E);
    n_mean = 2*pi/IP(i);
    omega_p(i) = -(3/2)*(R_E^2)*J2coeff(i)*n_mean*cosd(incl(i))/(a(i)*(1-ecc(i)^2))^2;
end
for q=1:length(r(1,:))
    rmag(q) = norm(r(:,q));
end
[rapo, apo_idx] = max(rmag);

plot3([0,r(1,apo_idx)],[0,r(2,apo_idx)],[0,r(3,apo_idx)],'m',LineWidth=2)
hold off

figure
plot(t,RAAN,'r')
grid on
title('Simulated Right Ascension of the Ascending Node vs Time')
xlabel('time')
ylabel('RAAN (\Omega) [deg]')

figure
plot(t,ecc,'b')
grid on
title('Simulated Eccentricity vs Time')
xlabel('time')
ylabel('e')

figure
plot(t,omega_p)
