function dcm_rtn = build_rtn_frame_fast(x_cart)
%BUILD_RTN_FRAME: Construct the RTN (radial, transverse, normal)
%frame associated with a given cartesian (inertial) state vector in the 2BP
%   Inputs:
%       x_cart:     Inertial (Cartesian) state vector
%   Outputs:
%       dcm_rtn:    nCr = [rhat thetahat hhat]
    r = x_cart(1:3);
    v = x_cart(4:6);
    h = cross(r, v);

    rhat = r/norm(r);
    hhat = h/norm(h);
    thetahat = cross(hhat, rhat);

    dcm_rtn = [rhat thetahat hhat];

end

