function inEclipse = eclipseCheck(rSun,rSat)
%=========================================================================%
% FUNCTION: eclipseCheck
% AUTHOR: Your Name (your email), 2021
% DESCRIPTION: Computes whether a satellite is in eclipse using the
% assumption of a spherical Earth and point-like Sun.
% INPUTS:
%   > rSun = [3x1] position vector for the Sun
%   > rSat = [3x1] position vector for the Satellite
% OUTPUTS:
%   > inEclipse = flag set to '1' if in eclipse and '0' otherwise.
%=========================================================================%

% Earth radius
Re = 6371; % [km]

% Check eclipse status
uhat = cross(rSun,rSat)/norm(cross(rSun,rSat));
phat = cross(uhat,(rSun - rSat)/norm(rSun - rSat));
if dot(rSat,phat) < Re && dot(rSat,rSun) < 0
    inEclipse = 1;
else
    inEclipse = 0;
end


end