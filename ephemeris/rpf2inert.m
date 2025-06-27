function [x_inert_nd] = rpf2inert(true_epochs, x_rpf, model, opts)
%RPF2INERT: Convert a set of non-dim RPF states to corresponding 
% NON-DIMENSIONAL states in the inertial J2000 frame.
%   Inputs:
%       true_epochs:    True epochs (seconds since J2000) assoc. with x_rpf
%       x_rpf:          Rotating-pulsating frame states (non-dim)
%       model:          EphemForceModel specifying conversion
%   Options:
%       uniform_time:   Flag to use uniform time (1) or non-uniform (0) in
%                       performing state conversion.  Default is uniform
%   Outputs:
%       x_inert:        Inertial (non-dim) states   
arguments
    true_epochs (1, :) double
    x_rpf (6, :) double
    model (1, 1) EphemForceModel
    opts.uniform_time (1, 1) {mustBeInteger} = 1;
end
    [M, bary_states] = rpf2inert_mat(true_epochs, model, "uniform_time",...
        opts.uniform_time);

    x_rpf = reshape(x_rpf, 6, 1, []);
    x_inert_dim = pagemtimes(M, x_rpf)+bary_states;

    x_inert_nd = model.nondimensionalize(x_inert_dim);
end

