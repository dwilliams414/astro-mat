%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Location of the L2 Lagrange Point
% Inputs:
% mu - gravitation parameter
% 
% Outputs:
% l2: column vector with location of L2
% 
% Validation Date: 09/27/2021 against Emily's Thesis with Lunar Mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l2 = findL2(mu, gamma20)
    f = @(gamma2) (1-mu)/(1+gamma2)^2+mu/(gamma2^2)-(1-mu+gamma2);
    df = @(gamma2) -2*(1-mu)/(1+gamma2)^3-2*mu/gamma2^3-1;
    
    [gamma2, numIter] = newtonRaphson1D(f, df, gamma20, 1e-12, 10000);
    l2 = (1-mu)+gamma2;
    l2 = [l2;0;0];
    % numIter
end