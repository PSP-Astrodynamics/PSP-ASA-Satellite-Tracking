function [a_m, e, inc_rad, RAAN_rad, wP_rad, nu] = getClassicalElements(tleLine2)

MU = 3.986004415e14; %Earth's gravitational parameter (m^3/s^2)
L2 = textscan(tleLine2,'%d%6d%f%f%f%f%f%f');

%% Orbital Elements
inc       = L2{1,3};                         % Inclination [deg]
inc_rad   = inc * pi / 180;                  % Inclination [rad]
RAAN      = L2{1,4};                         % Right Ascension of the Ascending Node [deg]
RAAN_rad  = RAAN * pi / 180;                 % Right ascension of the Ascending Node [rad]
e         = L2{1,5}/1e7;                     % Eccentricity 
wP        = L2{1,6};                         % Argument of periapsis [deg]
wP_rad    = wP * pi / 180;                   % Argument of periapsis [rad]
n         = L2{1,8};                         % Mean motion [Revs per day]
nu        = n * 7.2722052166431e-5;          % Mean motion [rad / sec]
a         = (MU/(n*2*pi/(24*3600))^2)^(1/3); % Semi-major axis [m]
a_m       = a / 1000;                        % Semi-major axis [km]

end