%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% monodromy3D: Computes the monodromy matrix based on a given half-period
% STM for x-axis symmetry
%
% Validated 11/3/2021
% Inputs:
% phiHalf: The half period STM
% 
% Outputs:
% M: The monodromy matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = monodromy3DxAxis(phiHalf)
    G = diag([1 -1 1 -1 1 -1]);
    Omega = [0 1 0; -1 0 0; 0 0 0];
    symp1 = [zeros(3) -eye(3); eye(3) -2*Omega];
    symp2 = [-2*Omega eye(3); -eye(3) zeros(3)];
    
    M = G*symp1*(phiHalf.')*symp2*G*phiHalf;
end