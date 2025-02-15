function vec_rotated = rotateXYZrth(vec,inc,raan,argp,nu,FROM,TO)
    
    theta = nu+argp;

    ICR = [cos(raan)*cos(theta) - sin(raan)*cos(inc)*sin(theta), -cos(raan)*sin(theta)...
    - sin(raan)*cos(inc)*cos(theta), sin(raan)*sin(inc);
       sin(raan)*cos(theta) + cos(raan)*cos(inc)*sin(theta),...
       -sin(raan)*sin(theta) + cos(raan)*cos(inc)*cos(theta), -cos(raan)*sin(inc);
       sin(inc)*sin(theta), sin(inc)*cos(theta), cos(inc)];
    
    if ((FROM=='RTH') & (TO=='XYZ'))
        vec_rotated = (ICR*vec')';
    elseif ((FROM=='XYZ') & (TO=='RTH'))
        vec_rotated = (ICR'*vec')';
    end
end

