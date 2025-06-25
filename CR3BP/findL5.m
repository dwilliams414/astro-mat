%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Location of the L5 Lagrange Point
% Inputs:
% mu - gravitation parameter
% 
% Outputs:
% l5: column vector with location of L5
% 
% Validation Date: 09/27/2021 against Emily's Thesis with Lunar Mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l5 = findL5(mu)
    l5 = [1/2-mu; -sqrt(3)/2; 0];
end