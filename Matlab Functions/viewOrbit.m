%=========================================================================%
% PROGRAM: viewOrbit.m
% AUTHOR: Justin Mansell (jmansell@purdue.edu), 2021
% DESCRIPTION: Creates a cool animation of the orbit for AAE 590.
% NOTE: You must run your solution for main.m
%=========================================================================%

% Create the globe
[fig, globe] = globe3D();

% Convert trajectory into ECEF
for ii = 1:numel(propJulianDates)
    [recef(:,ii),vecef(:,ii),~]  = convertECI2ECEF(propStates(1:3,ii),propStates(4:6,ii),propJulianDates(ii),Environment);
end

% Add to the globe
an = animatedline(recef(1,1)*1e3,recef(2,1)*1e3,recef(3,1)*1e3,'Color','r');
for ii = 1:numel(propJulianDates)
    addpoints(an,recef(1,ii)*1e3,recef(2,ii)*1e3,recef(3,ii)*1e3)
    drawnow
end