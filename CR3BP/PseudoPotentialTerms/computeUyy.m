function Uyy = computeUyy(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);
    
    Uyy = 1-(1-mu)/d^3-mu/r^3+3*(1-mu)*y^2/d^5+3*mu*y^2/r^5;
end