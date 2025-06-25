function E = meanAnomaly2eccentricAnomaly_rad(M, e, tol)
%MEANANOMALY2ECCENTRICANOMALY: Convert mean anomaly values (radians) to
%eccentric anomaly values (radians).
%   Inputs:
%       M:      Mean anomaly values (radians)
%       e:      Eccentricity values
%   Outputs:
%       E:      Eccentric anomaly values (radians)
    sz_init = size(M);
    M = M(:); e = e(:);
    E0 = M + (e-e.^3/8).*sin(M)+(e.^2/2).*sin(2*M)+(3*e.^3/8).*sin(3*M);

    E = zeros(size(M));
    for ii = 1:length(E)
        if (e(ii) >= 1)
            error("Invalid eccentricity!");
        end
        E(ii) = keplerSolveRaphson(M(ii), e(ii), E0(ii), tol, 100);
    end

    E = reshape(E, sz_init);
    
end

function [E, numIter, dE_i, dM_i] = keplerSolveRaphson(M, e, E0, dEtol, maxIter)
    numIter = 0;
    
    dM_i = M - (E0-e*sin(E0));
    dE_i = dM_i/(1-e*cos(E0));
    E = E0; % For the case when the initial guess is good enough
    
    while((abs(dE_i) > dEtol) && (numIter < maxIter))
        dM_i = M-(E-e*sin(E));
        dE_i = dM_i/(1-e*cos(E));
        E = E + dE_i;
        numIter = numIter + 1;
    end
    if numIter >= maxIter
        fprintf("Warning: Maximum Iterations Reached\n");
        fprintf("M: %.8f\n", M);
        fprintf("e: %.8f\n", e)
    end
end