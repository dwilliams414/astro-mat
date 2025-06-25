function Ux = computeUx(pos,mu)
%computeUx: Computes the partial of U* wrt x (CR3BP)
%   pos: nd Position vector (3x1)
%   mu: Mass paramter (nd)
    x = pos(1);
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);
    
    Ux = -(1-mu)*(x+mu)/d^3-mu*(x-1+mu)/r^3+x;
end