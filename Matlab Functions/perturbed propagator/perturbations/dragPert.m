function a_D_xyz = dragPert(x, y, z, vx, vy, vz, t, mu, S_m, Cd)
    
    R0 = 6378.137;  
    omega_E = 7.2921151467e-5; % rad/s

    r_vec = [x, y, z];
    r = norm(r_vec);

    v_vec = [vx, vy, vz];
    %v = norm(v_vec);

    alt = r-R0;

    [~, ~, rho] = airProperties(alt); % kg/km^3

    v_rel_vec = v_vec - cross(omega_E*[0,0,1], r_vec);
    %v_rel_vec = v_vec; % ignore earth rotation
    v_rel = norm(v_rel_vec);

    % Compute drag perturbation acceleration vector (km/s²)
    a_D_xyz = -0.5*rho*Cd*S_m*(v_rel^2)*v_rel_vec/norm(v_rel_vec);
    a_D_xyz = a_D_xyz';
end
