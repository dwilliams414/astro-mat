function uxxy = computeUxxy(pos, mu)
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uxxy = 3*(1-mu)*y/d^5+3*mu*y/r^5-15*(1-mu)*(x+mu)^2*y/d^7-15*mu*(x+mu-1)^2*y/r^7;
end