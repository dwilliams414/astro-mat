function F = F_apse(X, x0, tmin, tmax, sys)
%F_TA Summary of this function goes here

    tspan = [0 X(1)];
    
    
    [t, x] = IntegrateCR3BP(x0, tspan, sys.mu);

    F = zeros(3, 1);
    F(1) = dot(x(end, 1:3)', x(end, 4:6)');
    F(2) = (t(end)-tmin)-X(2)^2;
    F(3) = (tmax-t(end))-X(3)^2;
end

