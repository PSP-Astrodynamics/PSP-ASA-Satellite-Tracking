function [t, Xkep] = keplerianProp(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert)
    E0 = 2*atan(((1-ecc0)/(1+ecc0))^0.5 * tan(nu0/2));
    M0 = E0 - ecc0*sin(E0);
    
    Xkep0 = [a0,ecc0,inc0,raan0,argp0,M0];
    [t, Xkep] = ode45(@(t,Xkep) dXkepdt(t,Xkep,simT,mu,pert),simT,Xkep0,opt);
end

function dXkepdt = dXkepdt(t,Xkep,simT,mu,pert)
 
    a = Xkep(1);
    ecc = Xkep(2);
    inc = Xkep(3);
    raan = Xkep(4);
    argp = Xkep(5);
    M = Xkep(6);

    p = a*(1-ecc^2);
    n = sqrt(mu/a^3);
    h = sqrt(p*mu);
    b = a*sqrt(1-ecc^2);

    E = meanToEcc(M, ecc);
    nu = 2*atan(((1+ecc)/(1-ecc))^0.5 * tan(E/2));

    r = p/(1+ecc*cos(nu));

    f0 = [0,0,0,0,0,n]';

    Bx = (1/h)*[2*a^2*ecc*sin(nu), 2*a^2*p/r, 0;...
        p*sin(nu), (p+r)*cos(nu)+r*ecc, 0;...
        0, 0, r*cos(nu+argp);...
        0, 0, r*sin(nu+argp)/sin(inc);
        -p*cos(nu)/ecc, (p+r)*sin(nu)/ecc, -r*sin(nu+argp)/tan(inc);...
        b*p*cos(nu)/(a*ecc)-2*b*r/a, -b*(p+r)*sin(nu)/(a*ecc), 0];

    [x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,nu,mu);

    %a_J2_xyz = -(3*mu*J2*(R0^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    %(1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z];
    
    a_pert_xyz = [0; 0; 0];
    for jj=1:length(pert)
        a_pert_xyz = a_pert_xyz + pert{jj}(x,y,z,vx,vy,vz,t,mu);
        
    end
    %a_pert_xyz = double(a_pert_xyz);
    % a_pert_xyz1 = -(3*mu*1.0826e-3*(6378.1^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    % (1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z]

    a_pert_rth = rotateXYZrth(a_pert_xyz',inc,raan,argp,nu,'XYZ','RTH');
    a_pert_rth = a_pert_rth';

    dXkepdt = f0 + Bx*a_pert_rth;
end