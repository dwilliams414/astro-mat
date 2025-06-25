%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeJC: Computes the Jacobi Constant in the CR3BP
% Inputs:
% rho: position state [nd] (COLUMN vector)
% v: velocity state [nd] (COLUMN vector)
% Outputs:
% JC: Jacobi Constant
% Validated 09/27/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function JC = computeJC(rho, v, mu)
    d = sqrt((rho(1)+mu)^2+rho(2)^2+rho(3)^2);
    r = sqrt((rho(1)-1+mu)^2+rho(2)^2+rho(3)^2);
    UStar = (1-mu)/d+mu/r+0.5*(rho(1)^2+rho(2)^2);
    JC = 2*UStar - norm(v)^2;
end