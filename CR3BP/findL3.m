%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Location of the L3 Lagrange Point
% Inputs:
% mu - gravitation parameter
% 
% Outputs:
% l3: column vector with location of L3
% 
% Validation Date: 09/27/2021 against Emily's Thesis with Lunar Mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l3 = findL3(mu, gamma30)
    f = @(g) (1-mu)/g^2+mu/(1+g)^2-(mu+g);
    df = @(g) -2*(1-mu)/g^3-2*mu/(1+g)^3-1;
    
    [gamma3, numIter] = newtonRaphson1D(f, df, gamma30, 1e-12, 1000);
    l3 = -mu-gamma3;
    l3 = [l3;0;0];
    % numIter
end