function nMr6 = rotating2inertial_cr3bp_M(t)
%rotating2inertial_cr3bp: Transformation matrix relating a (Barycentric)
%CR3BP rotating frame state to (Barycentric) nondimensional inertial
%positioni and velocity.
%   Inputs:
%       t:              Time of evaluation (equivalent to theta in rad)
%   Outputs:
%5      nMr6
arguments
    t          (1, :) double
end
    nMr6 = zeros(6, 6, length(t));
    nCr = zeros(3, 3, length(t));

    % Build DCM - pagewise
    nCr(1, 1, :) = cos(t);
    nCr(1, 2, :) = -sin(t);
    nCr(2, 1, :) = sin(t);
    nCr(2, 2, :) = cos(t);
    nCr(3, 3, :) = ones(1, length(t));

    w_tilde = cross_matrix([0; 0; 1]); % Angular velocity vector

    nMr6(1:3, 1:3, :) = nCr;
    nMr6(4:6, 4:6, :) = nCr;
    nMr6(4:6, 1:3, :) = pagemtimes(nCr, w_tilde);

end