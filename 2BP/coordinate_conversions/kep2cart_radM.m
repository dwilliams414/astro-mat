function x_cart = kep2cart_radM(x_kep, GM, tol)
%KEP2CART_RADM: Convert Keplerian state x_kep = [a e i raan aop M] (with
%angles measured in radians) to a cartesian state
a = x_kep(1, :);
e = x_kep(2, :);
inc = x_kep(3, :);
raan = x_kep(4, :);
aop = x_kep(5, :);
M   = x_kep(6, :);

% Obtain true anomaly
ta = meanAnomaly2trueAnomaly_rad(M, e, tol);

% Obtain Rotation matrix to perifocal frame
dcm_perifocal = body2_313(raan, inc, aop);

% Obtain position magnitudes
p = a.*(1-e.^2);
r_mag = p./(1+e.*cos(ta));

% Obtain State Vectors: Postion and velocity (p58 BMW)
r_perifocal = r_mag.*[cos(ta); sin(ta); zeros(1, length(ta))];
v_perifocal = sqrt(GM./p).*[-sin(ta); (e+cos(ta)); zeros(1, length(ta))];

% Reshape for page
r_perifocal = reshape(r_perifocal, 3, 1, []);
v_perifocal = reshape(v_perifocal, 3, 1, []);

% Multiply by DCMs
r_inert = pagemtimes(dcm_perifocal, r_perifocal);
v_inert = pagemtimes(dcm_perifocal, v_perifocal);

r_inert = reshape(r_inert, 3, []);
v_inert = reshape(v_inert, 3, []);

x_cart = [r_inert; v_inert];
end

