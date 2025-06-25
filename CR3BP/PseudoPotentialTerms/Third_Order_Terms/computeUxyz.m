function uxyz = computeUxyz(pos, mu)
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    uxyz = -15*(1-mu)*(x+mu)*y*z/d^7-15*mu*(x+mu-1)*y*z/r^7;

end