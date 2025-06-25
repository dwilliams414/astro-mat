function x_cart = mil2cart(x_mil, GM, tol)
%MIL2CART: Convert the milankovitch state vector to cartesian coordinates
%   Inputs:
%       x_mil       [h; e; L] (radians)
%       GM          Mass parameter
%       tol:        Tolerance for solving for eccentric anomaly
%    Outputs:
%       x_cart:     Cartesian state vector
    x_kep = mil2kep_rad(x_mil, GM);
    x_cart = kep2cart_radM(x_kep, GM, tol);
end

