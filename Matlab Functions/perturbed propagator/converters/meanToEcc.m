function E = meanToEcc(M, e)
    tol = 1e-8;
    % Solve Kepler's equation for Ef
    E = M;  % Initial guess for E_f
    while abs(M - (E - e * sin(E))) > tol
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E));
    end
end