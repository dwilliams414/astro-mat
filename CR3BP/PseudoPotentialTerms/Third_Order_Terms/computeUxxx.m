function uxxx = computeUxxx(pos, mu)
    d = compute_d(pos, mu);
    r = compute_r(pos, mu);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    t1 = 6*(1-mu)*(x+mu)/d^5+6*mu*(x+mu-1)/r^5;
    t2 = 3*(1-mu)*(x+mu)/d^5+3*mu*(x+mu-1)/r^5;
    t3 = -15*(1-mu)*(x+mu)^3/d^7-15*mu*(x+mu-1)^3/r^7;

    uxxx = t1+t2+t3;
end