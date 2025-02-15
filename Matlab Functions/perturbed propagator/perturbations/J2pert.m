function a_J2_xyz = J2pert(x,y,z,vx,vy,vz,t,mu)
    R0 = 6378.1; % km
    J2 = 1.0826e-3; % ~

    r = sqrt(x^2+y^2+z^2);
    
    a_J2_xyz = -(3*mu*J2*(R0^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    (1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z];
end