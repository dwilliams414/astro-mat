function x_inert_nd = rotating2inertial_cr3bp_nd(t, x_cr3bp, opts)
%ROTATING2INERTIAL_CR3BP: Convert a CR3BP rotating frame state to inertial
%state representation.  Options specify shift of center for ineretial
%representation.
%inertial frame, etc.
%   Inputs:
%       t:      Time of evaluation (equivalent to angle)
%       x:      CR3BP nondimensional state vectors
%   Options:
%       'CenterND':     3x1 vector specifying center of inertial frame,
%                       with respect to system barycenter (nondim)
%   Outputs:
%       x_inert_nd:     Nondimensional position and velocity states
%                       expressed in the inertial frame.
arguments
    t   (1, :) double
    x_cr3bp (6, :) double
    opts.CenterND (3, 1) double = zeros(3, 1);
end
    nMr = rotating2inertial_cr3bp_M(t);

    x_cr3bp_shift = [x_cr3bp(1:3, :)-opts.CenterND; x_cr3bp(4:6, :)];
    x_cr3bp_shift = reshape(x_cr3bp_shift, 6, 1, []);

    x_inert_nd = pagemtimes(nMr, x_cr3bp_shift);
    x_inert_nd = reshape(x_inert_nd, 6, []);
end

