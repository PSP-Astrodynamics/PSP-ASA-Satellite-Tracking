function Xmk_struct = milankovitchProp(a0,ecc0,inc0,raan0,argp0,nu0,simT,mu,opt,pert)
    [hx0,hy0,hz0,eccx0,eccy0,eccz0,L0] = kep2milank(a0,ecc0,inc0,raan0,argp0,nu0,mu);
    h0_vec_xyz = [hx0,hy0,hz0];
    ecc0_vec_xyz = [eccx0, eccy0, eccz0];
    
    Xmk0 = [h0_vec_xyz, ecc0_vec_xyz, L0];
    
    Xmk_struct = ode45(@(t,Xmk) dXmkdt(t,Xmk,simT,mu,pert),simT,Xmk0,opt);
end

function dXmkdt = dXmkdt(t,X,T,mu,pert)
 
    h_vec_xyz = X(1:3);
    ecc_vec_xyz = X(4:6);
    L = X(7);

    h = norm(h_vec_xyz);
    ecc = norm(ecc_vec_xyz);

    p = (h^2)/mu;
    a = p/(1-ecc^2);

    [~,~,inc,raan,argp,nu] = milank2kep(h_vec_xyz,ecc_vec_xyz,L,mu);
    % hx = h_vec_xyz(1);
    % hy = h_vec_xyz(2);
    % hz = h_vec_xyz(3);
    % ex = ecc_vec_xyz(1);
    % ey = ecc_vec_xyz(2);
    % ez = ecc_vec_xyz(3);
    % Kep = MIL2Kep(hx, hy, hz, ex, ey, ez, L, mu);
    % inc = Kep(3);
    % raan = Kep(4);
    % argp = Kep(5);
    % nu = Kep(6);

    %[x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,nu,mu);
    [x,y,z,vx,vy,vz] = kep2cart(a,ecc,inc,raan,argp,nu,mu);

    r_vec_xyz = [x,y,z];
    r = norm(r_vec_xyz);

    v_vec_xyz = [vx,vy,vz];

    f0 = [0,0,0,0,0,0,h/r^2]';
    
    r_bar_xyz = crossMatrix(r_vec_xyz);
    v_bar_xyz = crossMatrix(v_vec_xyz);
    h_bar_xyz = crossMatrix(h_vec_xyz);

    Bx = [r_bar_xyz;...
        (1/mu)*(v_bar_xyz*r_bar_xyz-h_bar_xyz);...
        r_vec_xyz(3)/(h*(h+h_vec_xyz(3)))*h_vec_xyz']; 

    %a_pert_xyz = -(3*mu*J2*(R0^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    %(1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z];

    a_pert_xyz = [0; 0; 0];
    for jj=1:length(pert)
        a_pert_xyz = a_pert_xyz + pert{jj}(x,y,z,vx,vy,vz,t,mu);
    end

    % R0 = 6378.1; % km
    % J2 = 1.0826e-3; % ~
    % r = sqrt(x^2+y^2+z^2);
    % a_pert_xyz = -(3*mu*J2*(R0^2)/(2*r^5))*[(1-5*(z/r)^2)*x;...
    % (1-5*(z/r)^2)*y; (3-5*(z/r)^2)*z];
    dXmkdt = f0 + Bx*a_pert_xyz;
    
    % R0 = 6378.1; % km
    % J2 = 1.0826e-3; % ~
    % b_x2 = 1/mu * (v_bar_xyz * r_bar_xyz - h_bar_xyz);
    % b_x3 = dot([0,0,1]', r_vec_xyz) / (h * (h + dot([0,0,1]', h_vec_xyz))) * h_vec_xyz';
    % b_x = [r_bar_xyz; b_x2; b_x3];
    % J2_term = -3 * mu * J2 * R0^2/(2 * r^5);
    % Z = r_vec_xyz(3);
    % ad_N = J2_term .* [(1 - 5 * Z^2/r^2) * r_vec_xyz(1); (1 - 5 * Z^2/r^2) * r_vec_xyz(2); (3 - 5 * Z^2/r^2) * r_vec_xyz(3)]
    % dXmkdt = f0 +  b_x * ad_N;

    
end