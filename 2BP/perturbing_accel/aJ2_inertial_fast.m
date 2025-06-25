function aJ2 = aJ2_inertial_fast(x_inert, GM, J2, r0)
%AJ2_INERTIAL: Compute the acceleration due to J2, expressed in the
%inertial (cartesian) coordinate frame
%   Inputs:
%       x_inert:        Inertial cartesian state vector
%       GM:             Mass parameter of primary
%       J2:             J2 coefficient of primary
%       r0:             Equatorial radius of primary
%   Outputs:
%       aJ2:            Acceleration due to J2 perturbation
    x = x_inert(1);
    y = x_inert(2);
    z = x_inert(3);
    r = norm(x_inert(1:3));

    leading_coeff = -3*GM*J2*r0^2/(2*r^5);

    aJ2 = leading_coeff * [
                            (1-5*z^2/(r^2))*x;
                            (1-5*z^2/(r^2))*y;
                            (3-5*z^2/(r^2))*z
                          ];

end

