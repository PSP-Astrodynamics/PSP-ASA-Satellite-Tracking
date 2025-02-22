function Xcar_struct = cartesianProp(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert)
    [x0,y0,z0,vx0,vy0,vz0] = kep2cart(a0,ecc0,inc0,raan0,argp0,nu0,mu);
    Xcar0 = [x0,y0,z0,vx0,vy0,vz0];

    Xcar_struct = ode45(@(t,Xkep) dScdt(t,Xkep,simT,mu,pert),simT,Xcar0,opt);
end

function dXcardt = dScdt(t,Xcar,simT,mu,pert)
    v_xyz = Xcar(4:6);

    r_xyz = Xcar(1:3);
    r = norm(r_xyz);
    r_hat_xyz = r_xyz/r;

    x = r_xyz(1);
    y = r_xyz(2);
    z = r_xyz(3);
    vx = v_xyz(1);
    vy = v_xyz(2);
    vz = v_xyz(3);

    a_kepler_xyz = -mu*r_hat_xyz/r^2;

    % a_J2_xyz = -(3*mu*J2*(R0^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    %     (1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z];

    a_pert_xyz = [0; 0; 0];
    for jj=1:length(pert)
        a_pert_xyz = a_pert_xyz + pert{jj}(x,y,z,vx,vy,vz,t,mu);
    end

    dXcardt = [v_xyz; a_kepler_xyz+a_pert_xyz];
end