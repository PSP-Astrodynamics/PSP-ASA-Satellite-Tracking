function [a_SRP_xyz] = SRPpert(x,y,z,vx,vy,vz,t,mu)
    a_E = 149597898; % km
    A_m = 5.4e-6; % km^2/kg
    G0 = 1.02e14; % kgkm/s^2
    
    n_E = sqrt(mu/a_E^3);
    psi = n_E*t; % rad

    rE_vec_xyz = a_E*[cos(psi); sin(psi); 0];
    r_vec_xyz = [x;y;z];

    d_vec_xyz = rE_vec_xyz+r_vec_xyz;
    d = norm(d_vec_xyz);
    s_hat_xyz = d_vec_xyz/d;

    a_SRP_xyz = A_m*G0*s_hat_xyz/d^2;
end