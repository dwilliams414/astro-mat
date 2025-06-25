% Computes d in CR3BP equations.  Validated 10/18/2021
function d = compute_d(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    d = sqrt((x+mu)^2+y^2+z^2);
end