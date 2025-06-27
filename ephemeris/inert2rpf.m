function [x_rpf, true_epochs] = inert2rpf(t_since_epoch_nd, x_inert_nd, ...
    model, opts)
%inert2rpf: Convert non-dimensional J2000 states to corresponding RPF
%states.
%   Inputs:
%       t_since_epoch_nd:   ND time since epoch associated with model.
%                           Measured in terms of characteristic quantities
%                           associated with model
%       x_inert_nd:         ND J2000 states wrt P1
%       model:              EphemForceModel
%   Options:
%       uniform_time        Flag to use uniform time (1) or nonuniform (0).
%                           Default is uniform
%   Outputs:
%       x_rpf:              Rotating pulsating frame states
arguments
    t_since_epoch_nd (1, :) double
    x_inert_nd (6, :)
    model (1, 1) EphemForceModel
    opts.uniform_time (1, 1) {mustBeInteger} = 1;
end
    true_epochs = model.get_true_epochs(t_since_epoch_nd);

    [M, bary_states] = rpf2inert_mat(true_epochs, model, ...
        'uniform_time', opts.uniform_time);
    x_J2000_dim = model.dimensionalize(x_inert_nd);

    Xvec = reshape(x_J2000_dim-reshape(bary_states, 6, []), 6, 1, []);

    % x_rpf = reshape(pagemtimes(pageinv(M), Xvec), 6, []);
    for k = 1:size(M, 3)
        x_rpf(:, k) = inv(M(:, :, k))*Xvec(:, 1, k);
    end
end