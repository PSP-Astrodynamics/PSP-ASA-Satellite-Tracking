function [p,f,g,hf,k,L] = kep2modEqui(a,ecc,inc,raan,argp,nu)
    p = a*(1-ecc^2);
    f = ecc*cos(argp+raan);
    g = ecc*sin(argp+raan);
    hf = tan(inc/2)*cos(raan);
    k = tan(inc/2)*sin(raan);
    L = raan+argp+nu;
end

