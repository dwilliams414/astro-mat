%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Location of the L1 Lagrange Point
% Inputs:
% mu - gravitation parameter
% 
% Outputs:
% L1: column vector with location of L1
% 
% Validation Date: 09/27/2021 against Emily's Thesis with Lunar Mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l1 = findL1(mu, gamma10)
    f = @(gamma) -(1-mu)/(1-gamma)^2+mu/gamma^2+(1-mu-gamma);
    df = @(gamma) -2*(1-mu)/(1-gamma)^3-2*mu/gamma^3-1;
    
    [gamma1, numIter] = newtonRaphson1D(f, df, gamma10, 1e-12, 1000);
    l1 = 1-mu-gamma1;
    l1 = [l1;0;0];
    % numIter
end