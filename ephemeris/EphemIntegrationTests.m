% Simple Test for Ephemeris Integration
clc; clear all; close all;

%% Load Spice Kernels
if (ispc())
    tls_file = "../kernels/latest_leapseconds.tls.pc";
else
    tls_file = "../kernels/latest_leapseconds.tls";
end
spk_file = "../kernels/de440.bsp";
pck_file = "../kernels/gm_de440.tpc";

cspice_furnsh(char(tls_file));
cspice_furnsh(char(spk_file));
cspice_furnsh(char(pck_file));

%% Specify Epoch Time
% et = cspice_str2et('2025 APR 08 00:00:00.0000');
et = 636292869.1853861;

% Propagate
x0 = [0.6061276  0.68460448  0.23174715 -0.40289883  0.25669492 0.14728152]';
model = EphemForceModel('et0', et);

tic;
[traj, opts] = integrate_ephemeris_uniform(x0, [0 1], model, 'prop_stm', 0);
toc;

figure(1); hold on; grid on; axis equal;
plot_traj(traj);

%% Validate A Matrix
NNR = NumericalNewtonRaphson(@(x) opts.eom(0, x));
Anum = NNR.finite_diff(x0);
A_analytic = A_ephem(et, x0, model);

cspice_kclear;