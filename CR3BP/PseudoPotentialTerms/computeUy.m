function Uy = computeUy(pos, mu)
%COMPUTEUY Compute the partial of U* wrt y
%   pos: nd Position vector (3x1)
%   mu: mass parameter
    y = pos(2);
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    Uy = y-(1-mu)*y/d^3-mu*y/r^3;
end

