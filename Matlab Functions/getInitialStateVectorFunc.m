function [rxyz,velxyz,alt,JD] = getInitialStateVectorFunc(tleLine1,tleLine2)

%% Parameters
MU = 3.986004415e14;        %Earth's gravitational parameter (m^3/s^2)
RE = 6378.1363e3;           %Earth's equatorial radius (m)

%% Read TLE and get classical orbital elements

L1 = textscan(tleLine1,'%d%6d%*c%5d%*3c%*2f%f%f%5d%*c%*d%5d%*c%*d%d%5d',[1,9]);
yr = textscan(tleLine1,'%d%6d%*c%5d%*3c%*0f%f%f%5d%*c%*d%5d%*c%*d%d%5d');
yr = 2000+floor(yr{4}/1e3); %Year
doy = L1{4}; %Day of year
L2 = textscan(tleLine2,'%d%6d%f%f%f%f%f%f');

%Orbit elements
inc       = L2{1,3};                % Inclination [deg]
RAAN      = L2{1,4};                % Right Ascension of the Ascending Node [deg]
e         = L2{1,5}/1e7;            % Eccentricity 
wP        = L2{1,6};                % Argument of periapsis [deg]
MeanAnom  = L2{1,7};                % Mean anomaly [deg]
n         = L2{1,8};                % Mean motion [Revs per day]
a = (MU/(n*2*pi/(24*3600))^2)^(1/3);     % Semi-major axis [km]

% True anomaly
KepEq = @(x) deg2rad(MeanAnom) - x + e*sin(x); %Kepler's equation
options = optimset('Display','off','TolX',1e-10); %Options to suppress output and decrease tolerance
EccAnom = fzero(KepEq,MeanAnom*pi/180,options);        %Solve for the eccentric anomaly [rad]
TA = 2*atan2(tan(EccAnom/2)*sqrt(1+e),sqrt(1-e));   %True Anomaly [rad]
TA = TA*180/pi; % True anomaly [deg]


%% Convert classical oribtal elements to cartesian state
R = orbitDCM(inc,RAAN,TA+wP); %Create the direction cosine matrix
p = a*(1-e^2); %Semi latus rectum
h = sqrt(MU*p); %Orbital angular momentum
rOrb = [p/(1+e*cosd(TA));0;0]; %Position in orbital frame
velOrb = [MU*e/h*sind(TA);h/rOrb(1);0]; %Velocity in orbital frame
rxyz = R*rOrb/1e3;
velxyz = R*velOrb/1e3;

%% Print results

% Julian date
TLEepoch = datetime(yr,1,0,0,0,0) + days(doy);
JD = days(TLEepoch - datetime(2000,1,1,0,0,0)) + 2451544.5;
alt = norm(rxyz)-RE/1000;

%% Orbital DCM function
function R = orbitDCM(inc,omeg,th)
%Creates the direction cosine matrix for an orbit. The matrix can then
%be used to convert between orbit frame coordinates (rhat,theta_hat,h_hat)
%and inertial coordinates (x,y,z).
R = eye(3);
R(1,1) = cosd(omeg)*cosd(th)-sind(omeg)*cosd(inc)*sind(th);
R(2,1) = sind(omeg)*cosd(th)+cosd(omeg)*cosd(inc)*sind(th);
R(3,1) = sind(inc)*sind(th);
R(3,2) = sind(inc)*cosd(th);
R(3,3) = cosd(inc);
R(1,2) = -cosd(omeg)*sind(th)-sind(omeg)*cosd(inc)*cosd(th);
R(2,2) = -sind(omeg)*sind(th)+cosd(omeg)*cosd(inc)*cosd(th);
R(1,3) = sind(omeg)*sind(inc);
R(2,3) = -cosd(omeg)*sind(inc);
return
end

end