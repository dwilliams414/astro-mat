function JC = compute_JC_vec(states, sys_mu)
%COMPUTE_JC_VEC Vectorized computation of the Jacobi Constant
% Validation: 09/30/2022
    arguments
        states (:, 6) double
        sys_mu (1, 1) double
    end

    Ustar = compute_UStar_vec(states, sys_mu);
    vsqrd = states(:, 4).^2+states(:, 5).^2+states(:, 6).^2;

    JC = 2*Ustar-vsqrd;
end

