function x_mil = cart2mil_rad(x_cart, GM, tol)
%CART2MIL: Convert a cartesian state vector to the Milankovitch vector
%   Inputs:
%       x_cart:     Inertial (cartesian state vector)
%       GM:         Mass parameter
%       tol:        Tolerance for solving Kepler's equation for True
%                   anomaly from mean anomaly
%   Outputs:
%       x_mil:      [hvec; evec; L] (radians)
    x_kep = cart2kep_radM(x_cart, GM);
    ta_vals = meanAnomaly2trueAnomaly_rad(x_kep(end, :), x_kep(2, :), tol);

    L = x_kep(4, :) + x_kep(5, :) + ta_vals;

    % Angular Momentum vector
    hvec = cross(x_cart(1:3, :), x_cart(4:6, :), 1);

    % Eccentricity Vector
    evec = compute_eccentricity_vec(x_cart, GM);

    x_mil = [hvec; evec; mod(L, 2*pi)];
end

