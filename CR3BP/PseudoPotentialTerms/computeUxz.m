%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeUxz: Computes U*_xz for the CR3BP
% Validated 10/18/2021
% Inputs:
% pos: Position
% mu: Mass parameter
%
% Outputs:
% uxz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uxz = computeUxz(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);
    
    uxz = 3*(1-mu)*(x+mu)*z/d^5+3*mu*(x-1+mu)*z/r^5;
end