function [a,ecc,inc,raan,argp,nu] = milank2kep(h_vec_xyz,ecc_vec_xyz,L,mu)
    h = norm(h_vec_xyz);
    ecc = norm(ecc_vec_xyz);

    p = (h^2)/mu;
    a = p/(1-ecc^2);

    h_hat_xyz = h_vec_xyz/h;
    ecc_hat_xyz = ecc_vec_xyz/ecc;

    % 3D orb elements
    inc = acos(h_hat_xyz(3));

    % raan = [asin(h_hat_xyz(1)/sin(inc)), pi-asin(h_hat_xyz(1)/sin(inc)),...
    %     acos(-h_hat_xyz(2)/sin(inc)), 2*pi-acos(-h_hat_xyz(2)/sin(inc))];
    % raan = nonuniqueAngle(raan);
    
    z_hat_xyz = [0, 0, 1];
    nodes_hat_xyz = cross(z_hat_xyz, h_hat_xyz)/norm(cross(z_hat_xyz, h_hat_xyz));
    
    if (nodes_hat_xyz(2)>=0)
        raan = acos(nodes_hat_xyz(1));
    else
        raan = -acos(nodes_hat_xyz(1));
    end
    %raan = atan2(h_hat_xyz(1), -h_hat_xyz(2));
    
    if (ecc_hat_xyz(3)>=0)
        argp = acos(dot(nodes_hat_xyz, ecc_hat_xyz));
    else
        argp = -acos(dot(nodes_hat_xyz, ecc_hat_xyz));
    end
    %ecc_hat_xyzperp = cross(h_hat_xyz, ecc_hat_xyz)/norm(cross(h_hat_xyz, ecc_hat_xyz));
    %argp = atan2(ecc_hat_xyz(3), ecc_hat_xyzperp(3));

    nu = L-raan-argp;
end

