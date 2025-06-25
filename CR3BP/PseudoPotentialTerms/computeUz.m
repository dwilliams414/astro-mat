function Uz = computeUz(pos,mu)
%COMPUTEUZ Compute the partial of the pseudopotential U* wrt z
%   pos: ND position vector (3x1)
%   mu: Mass parameter of system
    z = pos(3);
    r = compute_r(pos, mu);
    d = compute_d(pos, mu);

    Uz = -(1-mu)*z/d^3-mu*z/r^3;
end

