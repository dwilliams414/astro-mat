function x_kep = equi2kep_radM(x_equi)
%EQUI2KEP_RAD: Convert modified equinoctial elements to keplerian elements, 
% assuming angles in radians and outputting mean anomaly for the 
% fast Keplerian
% variable
%   Inputs:
%       x_equi:     [p; f; g; h; k; L] (radians for L)
%   Outputs:
%       x_kep:      [a; e; inc; raan; aop; M] (radians for angles)

% Unpack equinoctial
    p = x_equi(1, :);
    f = x_equi(2, :);
    g = x_equi(3, :);
    h = x_equi(4, :);
    k = x_equi(5, :);
    L = x_equi(6, :);

    % Get Keplerian
    e = sqrt(f.^2+g.^2);
    a = p./(1-e.^2);
    inc = 2*atan(sqrt(h.^2+k.^2));
    raan = atan2(k, h);
    aop = atan2(g, f)-raan;
    ta = L-raan-aop;

    M = trueAnomaly2meanAnomaly_rad(ta, e);
    M = mod(M, 2*pi);

    x_kep = [a; e; inc; raan; aop; M];
    x_kep(4:end, :) = mod(x_kep(4:end, :), 2*pi);
end

