function Xmod_struct = modEqProp(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert)
    [p0,f0,g0,hf0,k0,L0] = kep2modEqui(a0,ecc0,inc0,raan0,argp0,nu0);
    
    Xmod0 = [p0,f0,g0,hf0,k0,L0];  
    
    Xmod_struct = ode45(@(t,Xmod) dXmoddt(t,Xmod,simT,mu,pert),simT,Xmod0,opt);
end

function dXmoddt = dXmoddt(t,s,T,mu,pert)
 
    p = s(1);
    f = s(2);
    g = s(3);
    h = s(4);
    k = s(5);
    L = s(6);

    q = 1+f*cos(L)+g*sin(L);
    s = sqrt(1+h^2+k^2);

    f0 = [0,0,0,0,0,sqrt(mu*p)*(q/p)^2]';

    Bx = sqrt(p/mu)*[0, 2*p/q, 0;...
        sin(L), ((q+1)*cos(L)+f)/q, -g*(h*sin(L)-k*cos(L))/q;...
        -cos(L), ((q+1)*sin(L)+g)/q, f*(h*sin(L)-k*cos(L))/q;...
        0, 0, (s^2/(2*q))*cos(L);...
        0, 0, (s^2/(2*q))*sin(L);...
        0, 0, (h*sin(L)-k*cos(L))/q];
    
    % bring back orbital elements
    a = p/(1-f^2-g^2);
    ecc = sqrt(f^2+g^2);
    inc = atan2(2*sqrt(h^2+k^2), 1-h^2-k^2);
    raan = atan2(k,h);
    argp = atan2(g*h-f*k, f*h+g*k);
    nu = L - atan(g/f);

    r = p/(1+ecc*cos(nu));

    [x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,nu,mu);

    % a_J2_xyz = -(3*mu*J2*(R0^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    % (1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z];

    a_pert_xyz = [0; 0; 0];
    for jj=1:length(pert)
        a_pert_xyz = a_pert_xyz + pert{jj}(x,y,z,vx,vy,vz,t,mu);
    end

    a_pert_rth = rotateXYZrth(a_pert_xyz',inc,raan,argp,nu,'XYZ','RTH');
    a_pert_rth = a_pert_rth';

    dXmoddt = f0 + Bx*a_pert_rth;
end