function traj = integrate_system(sys, x0, tspan, opts)
%INTEGRATE_SYSTEM: Integrator for a given dynamical system\
arguments
    sys (1, 1) function_handle
    x0 (:, 1) double
    tspan (1, :) double
    opts.abs_tol (1, 1) double = 1e-12;
    opts.rel_tol (1, 1) double = 1e-12;
    opts.integrator (1, 1) function_handle = @ode89;
end
    odeopts = odeset('AbsTol', opts.abs_tol, 'RelTol', opts.rel_tol);
    
    traj = opts.integrator(sys, tspan, x0, odeopts);
end

