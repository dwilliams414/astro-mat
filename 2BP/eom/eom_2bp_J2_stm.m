function dxdt = eom_2bp_J2_stm(t, x, GM, J2, r0)
%EOM_2BP_J2_STM: Equations of motion defining propagation of state
%transition matrix (STM) for spacecraft under Keplerian motion perturbed by
%J2
    dxdt = zeros(6+36, 1);

    % Perturbing acceleration
    a_pert = @(t, x, GM) aJ2_inertial_fast(x, GM, J2, r0);

    % State derivative
    dxdt(1:6) = eom2bp_cartesian(t, x, GM, a_pert);

    % STM
    phi = reshape(x(7:end), 6, 6);
    A = A_2bp_J2(t, x, GM, J2, r0);
    phidot = A*phi;

    dxdt(7:end) = phidot(:);
end

