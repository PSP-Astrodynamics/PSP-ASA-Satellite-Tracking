function [a,ecc,inc,raan,argp,nu] = modEqui2kep(p,f,g,hf,k,L)
    a = p/(1-f^2-g^2);
    ecc = sqrt(f^2+g^2);
    inc = atan2(2*sqrt(hf^2+k^2), 1-hf^2-k^2);
    raan = atan2(k,hf);
    argp = atan2(g*hf-f*k, f*hf+g*k);
    nu = L - atan(g/f);
end

