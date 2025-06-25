function uyyx = computeUyyx(pos, mu)
d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uyyx = 3*(1-mu)*(x+mu)/d^5+3*mu*(x+mu-1)/r^5-15*(1-mu)*(x+mu)*y^2/d^7-15*mu*(x+mu-1)*y^2/r^7;
end