function x_cr3bp = inertial2rotating_cr3bp_nd(t, x_inert_nd, ...
    opts)
%INERTIAL2ROTATING_CR3BP_ND Convert inertial states to the CR3BP sidereal
%rotating frame.  It is assumed that the states are appropriately
%nondimensionalized when input from the inertial frame, using the desired
%CR3BP parameters.
%   Inputs:
%       t:              Time values (non-dimensional!)
%       x_inert_nd:     Nondim inertial states
%   Options:
%       Center:         Vector specifying position offset to apply.  For
%                       example, if converting from P1 centered states,
%                       pass in the location of P1 in the rotating frame,
%                       i.e., [-mu; 0; 0]
%   Outputs:
%       x_cr3bp:        Standard CR3BP nondim analogs
arguments
    t       (1, :) double
    x_inert_nd (6, :) double
    opts.Center (3, 1) double = zeros(3, 1);
end
    x_inert_nd = reshape(x_inert_nd, 6, 1, []);

    M_mats = rotating2inertial_cr3bp_M(t);

    x_cr3bp = pagemtimes(pageinv(M_mats), x_inert_nd);
    x_cr3bp = reshape(x_cr3bp, 6, []);

    x_cr3bp(1:3, :) = x_cr3bp(1:3, :) + opts.Center;
end

