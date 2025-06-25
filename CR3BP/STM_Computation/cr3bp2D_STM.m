%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cr3bp2D_STM: Propagates STM along with 2D CR3BP equations of motion of
% use in targeting and other applications
% 
% Inputs:
% t: Time 
% state: 20x1 state vecotor such that:
%   state(1:4) = [x;y;xDot;yDot];
%   state(4:end) = phiUnwrapped (Obtained from Phi(:), which appends columns
%   one after another
% mu: mass paramter
% 
% Outputs
% stateDot: 20x1 state vector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stateDot = cr3bp2D_STM(t, state, mu)
    stateDot = zeros(20,1); % initialize for return
    posVel2D = state(1:4); % Position velocity state
    phiUnwrapped = state(5:20); % Phi state
    
    posVel3D = [posVel2D(1:2); 0; posVel2D(3:4); 0]; % get 3D state
    deriv3D = cr3bp_derivs(t, posVel3D, mu); % get 3D deriv
    
    stateDot(1:2) = deriv3D(1:2); % x and y
    stateDot(3:4) = deriv3D(4:5); % xDot and yDot
    
    A = A_cr3bp2D(t, posVel2D, mu);
    phi = reshape(phiUnwrapped, [4 4]);
    phiDot = A*phi;
    
    stateDot(5:20) = phiDot(:);
end