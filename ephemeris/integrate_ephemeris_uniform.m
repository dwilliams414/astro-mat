function [traj, opts] = integrate_ephemeris_uniform(x0, tspan, model, opts)
%INTEGRATE_EPHEMERIS: Build and integrate an ephemeris force model.
%Assumes kernels are already loaded.
%   Inputs:
%       x0:     Initial state (ndim, inertial)
%       tspan:  Time span for integration
%       model:  EphemForceModel to use.
%   Options:
%       abs_tol:    Absolute tolerance for integration (1e-12)
%       rel_tol:    Relative tolerance for integration (1e-12)
%       integrator: ode89 Default
%       lstar:      Characteristic length (if desired to override)
%       tstar:      Characteristic time (if desired to override)
%       prop_stm:   Flag to propagate STM
arguments
    x0 (6, 1) double
    tspan (1, :) double
    model (1, 1) EphemForceModel;
    opts.abs_tol = 1e-12;
    opts.rel_tol = 1e-12;
    opts.integrator = @ode89;
    opts.prop_stm = 0;
end

    % Configure Equations of Motion
    eom = @(t, y) f_ephem_uni_natural(t, y, model);

    % Propagate
    if opts.prop_stm
        phi0 = eye(6);
        x0 = [x0; phi0(:)];
    end
    traj = integrate_system(eom, x0, tspan, "abs_tol", opts.abs_tol, ...
        'rel_tol', opts.rel_tol, 'integrator', opts.integrator);

    opts.eom = eom;
end

