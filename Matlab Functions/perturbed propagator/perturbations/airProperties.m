function [P, T, rho] = airProperties(h) 

    % Note: h (kilometers)
    
    %% put in meters
    h = 1000*h;
    
    R = 287; % J/KgK
    g = 9.81; % m/s^2
    
    % START (sea level)
    P_0 = 1.013 * 10^5; % Pa
    T_0 = 288.15; % K
    
    % FIRST LAYER (At h = 11 km) 
    P_layer1 = 2.265 * 10^4; % Pa
    h_layer1 = 11 * 10^3; % m
    T_layer1 = 216.66; % K
    
    % SECOND LAYER (At h = 25 km)
    P_layer2 = 2.5273 * 10^3; % Pa
    h_layer2 = 25 * 10^3; % m
    T_layer2 = 216.66; % K
    
    % FINAL VALUE
    h_final = 30 * 10^3; % m
    T_final = 231.24; % K
    
    % Calculations
    a1 = (T_layer1 - T_0) / (h_layer1 - 0); % K/m
    a3 = (T_final - T_layer2) / (h_final - h_layer2); % K/m
    
    if (h < h_layer1)
        % GRADIENT (LESS THAN FIRST LAYER)
        T = T_0 + a1*(h - 0);
        P = P_0 * (T / T_0)^(-g / (a1*R));
        rho = P/ (R * T);
    
    elseif ((h >= h_layer1) && (h < h_layer2))
        % ISOTHERMAL (MORE THAN FIRST LAYER)
        T = T_layer1; % K
        P = P_layer1 * exp(-(g / (R*T)) * (h - h_layer1));
        rho = P / (R * T);
    
    elseif (h >= h_layer2)
        % GRADIENT (MORE THAN SECOND LAYER)
        T = T_layer2 + a3*(h - h_layer2);
        P = P_layer2 * (T / T_layer2)^(-g / (a3*R));
        rho = P/ (R * T);
    end

    %% put back in km
    P = P *(1000/1)^2; % N/m^2 to N/km^2
    % T in K
    rho = rho*(1000/1)^3; % kg/m^3 to kg/km^3

end