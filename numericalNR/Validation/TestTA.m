clc; clear all; close all;
em_sys = load("em_constants.mat");
orb_dat = load("L2SHalo_FullStandardized.mat");
load("PlanetaryConstants.mat", 'moon');
syn_month_nd = 29.53*3600*24/em_sys.t_star;
%%
IP = 2/9*syn_month_nd;
x0 = interp_perp_orb_IP(IP, orb_dat, em_sys, 10);

F = @(X) F_apse(X, x0, 7*IP/8, 9*IP/8, em_sys);
X0 = [8*IP/9, 0.25, 0.25];

NR = NumericalNewtonRaphson(F);
[X_sol, F_sol, iter_error] = NR.solve(X0);