function A2bp = A_2bp_natural(t, x, GM)
%A_2BP_NATURAL Jacobian of two-body equations of motion with no
%perturbation
%   Inputs:
%       t:      Time of evalutation
%       x:      State vector to evaluate
%       GM:     Mass parameter
%   Outputs:
%       A_2bp_natural: Natural jacobian for two-body problem
    A2bp = zeros(6, 6, class(x));
    A2bp(1:3, 4:end) = eye(3);

    rmag = norm(x(1:3));
    r = x(1:3);
    rhat = r/rmag;
    A2bp(4:6, 1:3) = -GM/rmag^3*(eye(3)-3*(rhat*rhat'));
end

