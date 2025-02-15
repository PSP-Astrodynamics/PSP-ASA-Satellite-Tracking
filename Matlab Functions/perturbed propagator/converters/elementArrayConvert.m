function Xout = elementArrayConvert(Xin,FROM,TO,mu)
       
    if ((FROM=='CAR') & (TO=='KEP'))
        x=Xin(:,1); y=Xin(:,2); z=Xin(:,3);
        vx=Xin(:,4); vy=Xin(:,5); vz=Xin(:,6); 
        
        a = zeros(length(Xin(:,1)),1);
        ecc = a; inc=a; raan=a; argp=a; nu=a;
        Xout = [a, a, a, a, a, a];
        for j=1:length(a)
            [a(j),ecc(j),inc(j),raan(j),argp(j),nu(j)] = cart2kep(x(j),y(j),z(j),vx(j),vy(j),vz(j),mu);
            Xout(j,1:6) = [a(j),ecc(j),inc(j),raan(j),argp(j),nu(j)];
        end
    elseif ((FROM=='KEP') & (TO=='CAR'))
        a=Xin(:,1); ecc=Xin(:,2); inc=Xin(:,3);
        raan=Xin(:,4); argp=Xin(:,5); nu=Xin(:,6); 
        
        x = zeros(length(Xin(:,1)),1);
        y = x; z=x; vx=x; vy=x; vz=x;
        Xout = [x, x, x, x, x, x];
        for j=1:length(x)
            [x(j),y(j),z(j),vx(j),vy(j),vz(j)] = kep2cart(a(j),ecc(j),inc(j),raan(j),argp(j),nu(j),mu);
            Xout(j,1:6) = [x(j),y(j),z(j),vx(j),vy(j),vz(j)];
        end
    elseif ((FROM=='MOD') & (TO=='KEP'))
        p=Xin(:,1); f=Xin(:,2); g=Xin(:,3);
        hf=Xin(:,4); k=Xin(:,5); L=Xin(:,6); 
        
        a = zeros(length(Xin(:,1)),1);
        ecc = a; inc=a; raan=a; argp=a; nu=a;
        Xout = [a, a, a, a, a, a];
        for j=1:length(a)
            [a(j),ecc(j),inc(j),raan(j),argp(j),nu(j)] = modEqui2kep(p(j),f(j),g(j),hf(j),k(j),L(j));
            Xout(j,1:6) = [a(j),ecc(j),inc(j),raan(j),argp(j),nu(j)];
        end
    elseif ((FROM=='MKE') & (TO=='KEP'))
        h_vec_xyz = Xin(:,1:3);
        ecc_vec_xyz = Xin(:,4:6);
        L = Xin(:,7);
        
        a = zeros(length(Xin(:,1)),1);
        ecc = a; inc=a; raan=a; argp=a; nu=a;
        Xout = [a, a, a, a, a, a];
        for j=1:length(a)
            [a(j),ecc(j),inc(j),raan(j),argp(j),nu(j)] = milank2kep(h_vec_xyz(j,:),ecc_vec_xyz(j,:),L(j),mu);
            Xout(j,1:6) = [a(j),ecc(j),inc(j),raan(j),argp(j),nu(j)];
        end

    end

end