function evecs = compute_eccentricity_vec(x_cart, GM)
%COMPUTE_ECCENTRICITY_VEC: Compute eccentricity vector from a 2BP cartesian
%state
%   Inputs:
%       x_cart:     Cartesian state vector (6xn)
%       GM:         Mass parameter (1xn or scalar)
%   Outputs:
%       evecs:      Eccentricity vector associated with all states
    rvecs = x_cart(1:3, :);
    vvecs = x_cart(4:6, :);

    r = vecnorm(rvecs, 2, 1);
    v = vecnorm(vvecs, 2, 1);
    evecs = 1./GM.*((v.^2-GM./r).*rvecs-dot(rvecs,vvecs, 1).*vvecs);
end

