function uxxz = computeUxxz(pos, mu)
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uxxz = 3*(1-mu)*z/d^5+3*mu*z/r^5-15*(1-mu)*(x+mu)^2*z/d^7-15*mu*(x+mu-1)^2*z/r^7;
end