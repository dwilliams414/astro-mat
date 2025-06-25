%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cr3bp_STM: Full 3D CR3BP propagation with STM
% 
% Validated: 11/03/2021
% Inputs
% t: Time of evaluation
% state: 42-element state vector
%   state(1:6) = [x;y;z;xDot;yDot;zDot];
%   state(7:42) = Unwrapped STM
% mu: Mass parameter
%
% Outputs:
% stateDot: 42 element derivative vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stateDot = cr3bp_STM(t, state, mu)
    posVel = state(1:6);
    phiUnwrapped = state(7:42);
    phi = reshape(phiUnwrapped, [6 6]);
    
    % Position velocity derivative
    posVelDeriv = cr3bp_derivs(t, posVel, mu);
    
    % Phi Derivative
    A = A_cr3bp(t, posVel, mu);
    phiDot = A*phi;
    
    stateDot = [posVelDeriv; phiDot(:)];   
end