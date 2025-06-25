%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeUStar: Caluclates the value of U* for a given position and mass
% parameter.
% 
% Inputs:
% pos: nd position
% mu: Mass paramter
% 
% Outputs
% U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = computeUStar(pos, mu)
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);
    x = pos(1);
    y = pos(2);
    
    U = (1-mu)/d+mu/r+0.5*(x^2+y^2);
end