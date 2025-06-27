% Conversion test for rpf2inert
clc; clear all; close all;
em_sys = load("em_constants.mat");
orb_dat = load("L2SHalo_FullStandardized.mat");
x0 = orb_dat.x0_array(:, 23);
IP = orb_dat.IP_vals(23);

time_conv = 0;
%% Kernels
cspice_furnsh('../kernels/latest_leapseconds.tls.pc');
cspice_furnsh('../kernels/de440.bsp');
cspice_furnsh('../kernels/gm_de440.tpc');
%% Propagate and Get CR3BP Rotating Frame States
[tcr3bp, xcr3bp] = IntegrateCR3BP(x0, [0 IP], em_sys.mu);

% Get Patch Points
t_dim = tcr3bp*em_sys.t_star;
t_patch = linspace(0, IP, 5)*em_sys.t_star;
xpatch = interp1(t_dim, xcr3bp, t_patch)';
%% Ephemeris
model = EphemForceModel();
x_inert_nd = rpf2inert(t_patch+model.et0, xpatch, model, "uniform_time", time_conv);

% Convert back as a check
xpatch_check = inert2rpf(t_patch/model.tstar, x_inert_nd, model, "uniform_time", ...
    time_conv);
%% Plot
figure(1); hold on; grid on; axis equal;
plot_traj_row(xcr3bp);
plot3(xpatch(1, :), xpatch(2, :), xpatch(3, :), 'ro');

%% Plot Ephemeris
t_nd_ephem = model.get_nd_time_since_epoch(t_patch+model.et0);

for k = 1:size(x_inert_nd, 2)-1
    x0k = x_inert_nd(:, k);
    tspank = [t_nd_ephem(k), t_nd_ephem(k+1)];

    traj = integrate_ephemeris_uniform(x0k, tspank, model);
    
    x_rpf_k = inert2rpf(traj.x, traj.y, model, "uniform_time", time_conv);
    plot_traj_row(x_rpf_k', 'Color', 'b');
end

%% Propagate and Add from Inertial Propagation

cspice_kclear;