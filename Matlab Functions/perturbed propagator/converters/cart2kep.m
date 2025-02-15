function [a,e,inc,raan,argp,f] = cart2kep(x,y,z,vx,vy,vz,mu)
% tested: works fine to 10^-11

    r_vec_xyz = [x,y,z];
    r = norm(r_vec_xyz);
    r_hat_xyz = r_vec_xyz/r;

    v_vec_xyz = [vx,vy,vz];
    v = norm(v_vec_xyz);

    h_vec_xyz = cross(r_vec_xyz,v_vec_xyz);
    h = norm(h_vec_xyz);
    h_hat_xyz = h_vec_xyz/h;

    theta_hat_xyz = cross(h_hat_xyz,r_hat_xyz);

    ICR = [r_hat_xyz', theta_hat_xyz', h_hat_xyz'];
    RCI = ICR';

    r_vec_rth = (RCI*r_vec_xyz')';

    v_vec_rth = (RCI*v_vec_xyz')';

    eps = 0.5*v^2 - mu/r;

    a = -mu/(2*eps);

    p = (h^2)/mu;

    e = sqrt(1-p/a);
    
    % outbound condition
    if (v_vec_rth(1) >= 0)
        f = acos((p/r-1)/e);
    else
        f = -acos((p/r-1)/e);
    end

    % 3D orb elements
    inc = acos(h_hat_xyz(3));

    raan = [asin(h_hat_xyz(1)/sin(inc)), pi-asin(h_hat_xyz(1)/sin(inc)),...
        acos(-h_hat_xyz(2)/sin(inc)), 2*pi-acos(-h_hat_xyz(2)/sin(inc))];
    raan = nonuniqueAngle(raan);

    theta = [asin(r_hat_xyz(3)/sin(inc)), pi-asin(r_hat_xyz(3)/sin(inc)),...
        acos(theta_hat_xyz(3)/sin(inc)), 2*pi-acos(theta_hat_xyz(3)/sin(inc))];
    theta = nonuniqueAngle(theta);

    argp = theta - f;

end

