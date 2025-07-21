% Ephemeris State Transition Tests
clc; clear all; close all;
uniTime = 0;

cspice_furnsh('../../data/de440.bsp');
cspice_furnsh('../../data/latest_leapseconds.tls.pc');
cspice_furnsh('../../data/gm_de440.tpc');

orbDat = load("L2SHalo_FullStandardized.mat");
load("PlanetaryConstants.mat", 'moon', 'earth');
massRatio = moon.mu/(moon.mu+earth.mu);
lstarCR = moon.a;
tStarCR = sqrt(moon.a^3/(moon.mu+earth.mu));

x0 = orbDat.x0_array(:, 23);
IP = orbDat.IP_vals(23);

% Specify Ephemeris Force Model
model = EphemForceModel;

% Propagate CR3BP State
propCR3BP = @(x0, tspan) integrate_system(@(t, x) cr3bp_derivs(t, x, massRatio), x0, tspan);

tspanCR = linspace(1, 1+IP, 101);
traj = propCR3BP(x0, tspanCR);
xCR3BP = deval(traj, tspanCR);

figure(1); hold on; grid on; axis equal;
plot3(traj.y(1, :), traj.y(2, :), traj.y(3, :));

% Get the Ephemeris Epochs
ephemEpochs = model.et0 + tspanCR * tStarCR; % seconds past J2K
ephemND = model.get_nd_time_since_epoch(ephemEpochs);

% Get Ephem Analog - Initial Condition
x0J2K = rpf2inert(ephemEpochs(1), x0, model, "uniform_time", uniTime);

% Integrate Ephem
propEphem = @(x0, tspan) integrate_system(@(t, x) f_ephem_uni_natural(t, x, model), ...
    x0, tspan);

trajEphem = propEphem(x0J2K, ephemND);

% Get the CR States
ephemStates = deval(trajEphem, ephemND);
ephemRotPuls = inert2rpf(ephemND, ephemStates, model, "uniform_time", uniTime);
plot3(ephemRotPuls(1, :), ephemRotPuls(2, :), ephemRotPuls(3, :));


trajEphem.y(:, end)