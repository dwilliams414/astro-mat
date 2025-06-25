% Computes r in CR3BP equations.  Validated 10/18/2021
function r = compute_r(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    r = sqrt((x-1+mu)^2+y^2+z^2);
end