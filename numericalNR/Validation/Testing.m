clc; clear all; close all;

em_sys = load("em_constants.mat");
orb_dat = load("L2SHalo_FullStandardized.mat");
orb_i = 101;

X0 = orb_dat.x0_array(:, orb_i);

NR = NumericalNewtonRaphson(@(y) cr3bp_derivs(0.0, y, em_sys.mu));
DF = NR.finite_diff(X0);

DF_analytic = A_cr3bp(0, X0, em_sys.mu);

%% Test solving Simple equation
f = @(X) [X(1).^2-1; X(2).^3+1];

NR2 = NumericalNewtonRaphson(f, error_on_fail=1, err_tol=1e-14);
[X_sol, F_sol, iter_error] = NR2.solve([1.2; 2.2]);