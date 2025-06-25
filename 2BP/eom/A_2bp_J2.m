function A = A_2bp_J2(t, x, GM, J2, r0)
%A_2BP_J2: Jacobian matrix for two-body problem, perturbed by a J2
%accelration.
%   Inputs:
%       t:      Time of evaluation
%       x:      State vector (6x1) to evaluate
%       GM:     Two-body mass parameter
%       J2:     J2 coefficient
%       r0:     Radius of primary
    A_natural = A_2bp_natural(t, x, GM);
    aJ2_grad = gradJ2_accel(x(1:3), r0, GM, J2);

    A = A_natural + [zeros(3); eye(3)]*[aJ2_grad zeros(3)];
end

