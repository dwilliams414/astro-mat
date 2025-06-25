%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A_cr3bp2D: 2D Jacobian (Analytic) for CR3BP equations of motion in plane
%
% Inputs;
% t: Time to evaluate
% state: Position-velocity State (4x1) 
% mu: Mass parameter
% 
% Outputs:
% A: Jacobian Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = A_cr3bp2D(t, state, mu)
    pos = [state(1:2); 0];
    uxx = computeUxx(pos, mu);
    uyy = computeUyy(pos, mu);
    uxy = computeUxy(pos, mu);
    
    U = [uxx, uxy; uxy uyy];
    Omega = [0 2; -2 0];
    
    A = [zeros(2), eye(2); U, Omega];
    
end