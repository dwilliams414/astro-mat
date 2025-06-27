function traj = integrate_system(sys, x0, tspan, opts)
%INTEGRATE_SYSTEM: Integrate a given dynamical system, using MATLAB's
%built-in integrators
%   Inputs:
%       sys:    @(t, y) Function describing ODEs
%       x0:     Initial state vector
%       tspan:  Timespan for integration
%   Options:
%       'abs_tol':      Absolute tolerance for integration (1e-12 Default)
%       'rel_tol':      Relative tolerance for integration (1e-12 Default)
%       'integrator':   Handle for MATLAB propagator (@ode89 Default)
%   Outputs:
%       traj:           MATLAB ODE trajectory object
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

