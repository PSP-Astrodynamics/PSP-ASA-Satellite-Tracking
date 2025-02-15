function [a_moon_xyz] = moonPert(x,y,z,vx,vy,vz,t,GM_bary)
    mu_M = 4.9028e3; % km^3/s^2
    d_EM = 3.8475e5; % km
    
    n_M = sqrt(GM_bary/d_EM^3);

    psi = n_M*t; % rad

    rM_vec_xyz = d_EM*[cos(psi); sin(psi); 0];
    r_vec_xyz = [x;y;z];

    r_MS_vec_xyz = r_vec_xyz-rM_vec_xyz;
    a_moon_xyz = -mu_M*(r_MS_vec_xyz/norm(r_MS_vec_xyz)^3+...
        rM_vec_xyz/norm(rM_vec_xyz)^3);
end