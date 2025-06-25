function uyyy = computeUyyy(pos, mu)
d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uyyy = 9*(1-mu)*y/d^5+9*mu*y/r^5-15*(1-mu)*y^3/d^7-15*mu*y^3/r^7;