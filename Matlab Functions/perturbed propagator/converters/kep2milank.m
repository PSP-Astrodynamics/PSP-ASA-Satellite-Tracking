function [hx,hy,hz,eccx,eccy,eccz,L] = kep2milank(a,ecc,inc,raan,argp,nu,mu)

    [x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,nu,mu) ;
    r_vec_xyz = [x,y,z];
    v_vec_xyz = [vx,vy,vz];

    h_vec_xyz = cross(r_vec_xyz,v_vec_xyz);
    hx = h_vec_xyz(1);
    hy = h_vec_xyz(2);
    hz = h_vec_xyz(3);

    ecc_vec_xyz = cross(v_vec_xyz,h_vec_xyz)/mu - r_vec_xyz/norm(r_vec_xyz);
    eccx = ecc_vec_xyz(1);
    eccy = ecc_vec_xyz(2);
    eccz = ecc_vec_xyz(3);
    
    L = raan+argp+nu;
end

