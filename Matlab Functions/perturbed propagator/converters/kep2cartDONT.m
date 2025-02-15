 function [x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,nu,mu) 
    % tested: works fine to 10^-12

    theta = nu + argp;

    p = a * (1-ecc^2);
    h = sqrt(mu*p);

    r_vec_rth = p/(1+ecc*cos(theta - argp)) * [1, 0, 0];

    ICR = [cos(raan)*cos(theta) - sin(raan)*cos(inc)*sin(theta), -cos(raan)*sin(theta)...
    - sin(raan)*cos(inc)*cos(theta), sin(raan)*sin(inc);
       sin(raan)*cos(theta) + cos(raan)*cos(inc)*sin(theta),...
       -sin(raan)*sin(theta) + cos(raan)*cos(inc)*cos(theta), -cos(raan)*sin(inc);
       sin(inc)*sin(theta), sin(inc)*cos(theta), cos(inc)];

    r_vec_xyz = (ICR*r_vec_rth');
    r = norm(r_vec_xyz);

    x = r_vec_xyz(1);
    y = r_vec_xyz(2);
    z = r_vec_xyz(3);
    
    fdot = h/r^2;

    if (nu <= pi)
        rdot = ecc*sqrt(mu/p)*sin(nu);
    else
        rdot = -ecc*sqrt(mu/p)*sin(nu);
    end

    v_vec_rth = [rdot,r*fdot,0];
    v_vec_xyz = (ICR*v_vec_rth')';

    vx = v_vec_xyz(1);
    vy = v_vec_xyz(2);
    vz = v_vec_xyz(3);
end