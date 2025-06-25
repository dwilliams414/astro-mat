function ad_rtn = a_inert2rtn(t, x_cart, GM, a_fcn)
%A_INERT2RTN: Convert an acceleration vector from the inertial frame to the
%RTN frame.
%   Inputs:
%       t: Time of evaluation
%       x_cart: Cartesian state representation
%       GM:     Mass parameter of system
%       a_fcn:  a(t, x_cart, GM) -> acceleration in inertial frame
%   Outputs:
%       ad_rtn: Accleration in RTN frame
    a_inert = a_fcn(t, x_cart, GM);
    nCr = build_rtn_frame_fast(x_cart);

    ad_rtn = nCr'*a_inert;
end

