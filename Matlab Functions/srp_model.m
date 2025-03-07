clear; clc; close all

%% Initializations
a0 = 6.8155e3; %semi major axis [km]
ecc0 = 7.148e-4; %eccentricity
inc0 = deg2rad(97.3566); %inclination [rad]
raan0 = deg2rad(95.9229); %Right Ascension of the Ascending Node [rad]
argp0 = deg2rad(342.7667); %Argument of periapsis [rad]
nu0 = 2*pi*15.4297/86400; %mean motion [rad / sec]
simT = 0:0.01:86400*7; %simulation time [s]
mu = 3.986004415e5; %Earth gravitational parameter [km^3/s^2]
opt = odeset('RelTol',1e-13, 'AbsTol', 1e-13); %set options for ODE solver
pert = {@noPert,@SRPpert};
R = 6378; %Earth radius for plot [km]

%% Calculations/Function Calls
output = cartesianProp(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert);

t=output.x;
X=output.y;

figure(1)
%plotStaticPropagation(X')
%plot3(X(1,:),X(2,:),X(3,:))
plot3(X(1,:),X(2,:),X(3,:), 'r')
hold on
earthy(R)
axis equal

