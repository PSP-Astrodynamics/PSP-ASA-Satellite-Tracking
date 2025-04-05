%=========================================================================%
% PROGRAM: srp_model.m
% AUTHOR: Evan Paul, Krish Mehta, 2025
% DESCRIPTION: This file provides orbit propagation and related info 
% for the Satellite Tracking Project. 
%=========================================================================%
clear; clc; close all

%--------------------------USER SETTINGS----------------------------------%

%% Define Satellite TLE - UPDATE THESE VALUES
tleLine1 = '1 51085U 22002DF  25032.54719027  .00048242  00000+0  10702-2 0  9993';
tleLine2 = '2 51085  97.3538 102.9956 0006702 315.2737  44.7970 15.43560689169443';

%% Define Satellite Properties - UPDATE THESE VALUES
mass = 2;      % Satellite mass [kg]
Aref = .03405 / 1000000; % Satellite reference area [km^2]
Cd   = 0.2;    % Satellite drag coefficient [unitless]
S_m  = Aref / mass; % Satellite ballistic coefficient

%% Initial state vector values
[rxyz,velxyz,alt,JD1] = getInitialStateVectorFunc(tleLine1,tleLine2);

%% General - UPDATE THESE VALUES
simT = 0:10:86400*300; %simulation time [s]

%------------------------END USER SETTINGS--------------------------------%

%% Initializations

[a0,ecc0,inc0,raan0,argp0,nu0] = getClassicalElements(tleLine2);

mu = 3.986004415e5; % Earth gravitational parameter [km^3/s^2]
opt = odeset('RelTol',1e-13, 'AbsTol', 1e-13); % Set options for ODE solver
pert = {@noPert,@SRPpert,@J2pert,@dragPert};

%% Calculations/Function Calls
output = cartesianProp(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert,S_m,Cd,JD1);

t=output.x;
X=output.y;

figure(1)
plotStaticPropagation(X')
