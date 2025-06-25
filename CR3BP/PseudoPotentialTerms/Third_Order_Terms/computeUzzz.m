function uzzz = computeUzzz(pos, mu)
d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uzzz = 9*(1-mu)*z/d^5+9*mu*z/r^5-15*(1-mu)*z^3/d^7-15*mu*z^3/r^7;
end
