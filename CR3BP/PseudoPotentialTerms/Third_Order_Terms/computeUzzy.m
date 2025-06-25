function uzzy = computeUzzy(pos, mu)
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uzzy = 3*(1-mu)*y/d^5+3*mu*y/r^5-15*(1-mu)*y*z^2/d^7-15*mu*y*z^2/r^7;
end