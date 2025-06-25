%% Validation File For Third Order Partials
clc; clear all; close all;
load("em_constants.mat");
pos = 3*rand(3,1)-1.5;

% Validate Numeric Differencing against what we know works
uxxx_analytic = computeUxxx(pos, mu);
uxxx_num = fd(pos, @computeUxx, 1, mu);

uxxy_num = fd(pos, @computeUxx, 2, mu);
uxxy_analytic = computeUxxy(pos, mu);

uxxz_analytic = computeUxxz(pos, mu);
uxxz_num = fd(pos, @computeUxx, 3, mu);

uxyz_analytic = computeUxyz(pos, mu);
uxyz_num = fd(pos, @computeUxy, 3, mu);

uyyy_analytic = computeUyyy(pos, mu);
uyyy_num = fd(pos, @computeUyy, 2, mu);

uyyx_num = fd(pos, @computeUyy, 1, mu);
uyyx_analytic = computeUyyx(pos, mu);

uyyz = computeUyyz(pos, mu);
uyyz_num = fd(pos, @computeUyy, 3, mu);

uzzz = computeUzzz(pos, mu)
uzzz_num = fd(pos, @computeUzz, 3, mu)

uzzx = computeUzzx(pos, mu);
uzzx_num = fd(pos, @computeUzz, 1, mu);

uzzy = computeUzzy(pos, mu)
uzzy_num = fd(pos, @computeUzz, 2, mu)

function numeric = fd(pos, func, pos_index, mu)
    x0_plus = pos;
    x0_minus = pos;

    h = sqrt(eps);

    x0_plus(pos_index) = pos(pos_index)+h;
    x0_minus(pos_index) = pos(pos_index)-h;

    fplus = func(x0_plus, mu);
    fminus = func(x0_minus, mu);

    numeric = (fplus-fminus)/(2*h);
end