% Ephemeris Integration in Matlab
clc; clear all; close all;

% Furnish Kernels to Build Model
kernel_dir = ['..', filesep, 'data', filesep];
bsp_file = 'de440.bsp';
lsk_file = 'latest_leapseconds.tls.pc';
pck_file = 'gm_de440.tpc';

cspice_furnsh([kernel_dir, bsp_file]);
cspice_furnsh([kernel_dir, lsk_file]);
cspice_furnsh([kernel_dir, pck_file]);

% Build Ephem Model
model = EphemForceModel;

% Specify Propagation function builder
prop_builder = @() build_ephem_prop(model, kernel_dir, bsp_file, lsk_file,...
    pck_file);

prop_constant = parallel.pool.Constant(prop_builder);

% Specify Initial Condition Array
x0_array = randn(6, 1001);
tspan = linspace(0, 2*pi, 201);
xf_array = zeros(size(x0_array));
parfor ii = 1:size(x0_array, 2)
    ii
    traj = prop_constant.Value(x0_array(:, ii), tspan);
    xf_array(:, ii) = traj.y(:, end)';
end


%% Functions
function prop = build_ephem_prop(model, kernel_dir, bsp_name, lsk_name,...
    pck_name)
    
    % Load Kernels if not already loaded
    if cspice_ktotal('all') == 0
        cspice_furnsh([kernel_dir, bsp_name]);
        cspice_furnsh([kernel_dir, lsk_name]);
        cspice_furnsh([kernel_dir, pck_name]);
    end

    eom = @(t, y) f_ephem_uni_natural(t, y, model);
    

    prop = @(x0, tspan) integrate_system(eom, x0, tspan);
end