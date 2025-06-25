function x_kep = mil2kep_rad(x_mil, GM)
%MIL2KEP: Convert the Milankovitch state vector to the Keplerian vector
%   Inputs:
%       x_mil:      [h; e; L] (radians for L)
%       GM:         Mass parameter
%   Outputs:
%       x_kep:      [a; e; i; raan; aop; M] (radians)
    h_vec = x_mil(1:3, :);
    e_vec = x_mil(4:6, :);
    L = x_mil(end, :); % Radians

    h = vecnorm(h_vec, 2, 1);
    hhat = h_vec./h;

    % Compute eccentricity
    e = vecnorm(e_vec, 2, 1);
    ehat = e_vec./e;

    % Compute sma
    sma = h.^2./(GM.*(1-e.^2));

    % Compute Inclination
    inc = acos(hhat(3, :));

    % Compute Line of Nodes
    Zhat = repmat([0; 0; 1], 1, size(x_mil, 2));
    n_vec = cross(Zhat, h_vec, 1);
    nhat = n_vec./vecnorm(n_vec, 2, 1);

    % Get the RAAN
    raan = acos(nhat(1, :));
    shift_indices = nhat(2, :) < 0;
    raan(shift_indices) = 2*pi-raan(shift_indices);

    % Get the AOP
    aop = acos(dot(ehat, nhat, 1));
    shift_indices = ehat(3, :) < 0;
    aop(shift_indices) = 2*pi-aop(shift_indices);

    % Get the true anomaly
    ta = L-aop-raan;
    ta = mod(ta, 2*pi);

    % Convert to mean anomaly
    M = trueAnomaly2meanAnomaly_rad(ta, e);

    x_kep = [sma; e; inc; raan; aop; M];
end

