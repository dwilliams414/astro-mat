function dxdt =  eom2bp_cartesian(t, x, GM, ad)
%EOM2BP_CARTESIAN: EOM for integration of the 2BP equations of motion with
%a disturbing acceleration, ad
%   Inputs:
%       t:      Time of evaluation
%       x:      State vector cartesian
%       GM:     Mass parameter
%       ad:     ad(t, x, GM) defining disturbing acceleration function
%   Outputs:
%       dxdt:   Cartesian state vector derivative
    % for debugging
    dxdt = zeros(6, 1, class(x));
    % dxdt = zeros(6, 1);

    r_vec = x(1:3);
    r     = norm(r_vec);
    dxdt(1:3) = x(4:6);
    dxdt(4:end) = -GM*r_vec./(r^3)+ad(t, x, GM);
end

