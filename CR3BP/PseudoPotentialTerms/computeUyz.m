%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeUyz: Computes U*_yz for the CR3BP
% Validated 10/18/2021
% Inputs:
% pos: Position
% mu: Mass parameter
%
% Outputs:
% uyz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uyz = computeUyz(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);
    
    uyz = 3*(1-mu)*y*z/d^5+3*mu*y*z/r^5;
end