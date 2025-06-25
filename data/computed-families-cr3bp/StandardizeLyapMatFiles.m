% Creates a better saved file for the L1/L2 Lyapunov Orbits (no more
% annoying .r0, .v0 stuff)
clc; clear all; close all;

%% Load System Prameters
em_sys = load("em_constants.mat");
load("L1Lyap.mat");

figure(1); hold on; grid on;
P = Propagator(@cr3bp_derivs, em_sys.mu);
for i = 1:length(IP_vals)
    IC_vals{i} = [IC_vals{i}.r0; IC_vals{i}.v0];
    x0_array(:, i) = IC_vals{i};

    sol = P.propagate(x0_array(:, i), [0, IP_vals(i)]);
    plot_traj(sol);
end

clear sol P i em_sys
save("L1LyapStandardized.mat");