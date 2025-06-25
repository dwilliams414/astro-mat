function aJ2_rtn = aJ2_rtn_fast(x_cart, GM, J2, r0)
%AJ2_RTN_FAST: Obtain the J2 acceleration expressed in terms of the
%rotating RTN frame
%   Inputs:
%       x_cart:     Inertial (cartesian) state vector
%       GM:         Primary mass parameter
%       J2:         J2 coefficient
%       r0:         Equatorial radius of primary
%   Outputs:
%       aJ2_rtn:    J2 acceleration expressed in the RTN frame
    aJ2_inert = aJ2_inertial_fast(x_cart, GM, J2, r0);
    nCr = build_rtn_frame_fast(x_cart);

    aJ2_rtn = nCr'*aJ2_inert;
end

