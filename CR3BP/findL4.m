%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Location of the L4 Lagrange Point
% Inputs:
% mu - gravitation parameter
% 
% Outputs:
% l4: column vector with location of L1
% 
% Validation Date: 09/27/2021 against Emily's Thesis with Lunar Mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l4 = findL4(mu)
    l4 = [0.5-mu; sqrt(3)/2;0];
end