function uxy = computeUxy(pos, mu)
    x = pos(1);
    y = pos(2);
    z = pos(3);
    
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);
    
    uxy = 3*(1-mu)*(x+mu)*y/d^5 + 3*mu*(x-1+mu)*y/r^5;
end