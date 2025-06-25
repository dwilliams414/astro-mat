%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeUzz: Computes U*_zz for the CR3BP
% Validated 10/18/2021
% Inputs:
% pos: Position
% mu: Mass parameter
%
% Outputs:
% uzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uzz = computeUzz(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    r = compute_r(pos, mu);
    d = compute_d(pos, mu);
    
    uzz = -(1-mu)/d^3-mu/r^3+3*(1-mu)*z^2/d^5+3*mu*z^2/r^5;
end