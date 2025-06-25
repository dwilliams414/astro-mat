function x_equi = kep2equi_radM(x_kep, tol)
%KEP2EQUI_RAD: Convert a Keplerian state vector to modified equinoctial
%elements
%   Inputs:
%       x_kep:      [a; e; i; raan; aop; M] (radians for angles)
%   Outputs:
%       x_equi:     [p; f; g; h; k; L] % L in radians
    a = x_kep(1, :);
    e = x_kep(2, :);
    inc = x_kep(3, :);
    raan = x_kep(4, :);
    aop  = x_kep(5, :);
    M    = x_kep(6, :); % Mean anomaly in radians
    ta   = meanAnomaly2trueAnomaly_rad(M, e, tol);

    p = a.*(1-e.^2);
    f = e.*cos(aop + raan);
    g = e.*sin(aop + raan);
    h = tan(inc/2).*cos(raan);
    k = tan(inc/2).*sin(raan);

    L = mod(ta+raan+aop, 2*pi);

    x_equi = [p; f; g; h; k; L];
    
end

