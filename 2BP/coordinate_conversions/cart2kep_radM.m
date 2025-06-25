function x_kep = cart2kep_radM(x_cart, GM)
%CART2KEP: Convert Cartesian state vectors to the Keplerian equivalent.
%Note that the angles in the Keplerian state vector are reported in
%RADIANS and the fast variable is mean anomaly (radians)
%speed, so don't break it!
%   Inputs:
%       x_cart:     Cartesian state vectors (3xn)
%       GM:         Mass parameters (1xn) assumed (or scalar)
%   Outputs:
%       x_kep:      Keplerian state vector xk = [a e i raan aop M]

% Get angular momentum vector
r_vec   = x_cart(1:3, :);
v_vec   = x_cart(4:end, :);
h_vec   = cross(r_vec, v_vec, 1);

% Magnitudes
r       = vecnorm(r_vec, 2, 1);
v       = vecnorm(v_vec, 2, 1);
h       = vecnorm(h_vec, 2, 1);

% Unit Vectors
rhat = r_vec./r;
vhat = v_vec./v;
hhat = h_vec./h;

% Compute Eccentricity Vector (pg. 50 BMW) & value
% e_vec = 1./GM.*((v.^2-GM./r).*r_vec - dot(r_vec, v_vec, 1).*v_vec);
e_vec = compute_eccentricity_vec(x_cart, GM);
e = vecnorm(e_vec, 2, 1);
ehat = e_vec./e;

% Compute Line of Nodes
Zhat = repmat([0; 0; 1], 1, size(x_cart, 2));
n_vec = cross(Zhat, h_vec, 1);
nhat = n_vec./vecnorm(n_vec, 2, 1);

% Compute eccentricity
e     = vecnorm(e_vec, 2, 1);

% Compute semi-major axis
p     = h.^2./GM;
sma   = p./(1-e.^2);

% Compute Inclination
inc   = acos(hhat(3, :));

% Compute Right Ascension
Xhat = repmat([1; 0; 0], 1, size(x_cart, 2));
raan = acos(nhat(1, :));
shift_indices = nhat(2, :) < 0;
raan(shift_indices) = 2*pi-raan(shift_indices);

% Compute AOP
aop = acos(dot(nhat, ehat));
shift_indices = ehat(3, :) < 0;
aop(shift_indices) = 2*pi-aop(shift_indices);

% Compute True Anomaly
ta = acos(dot(rhat, ehat, 1));
shift_indices = dot(r_vec, v_vec, 1) < 0;
ta(shift_indices) = 2*pi-ta(shift_indices);

% Compute Mean Anomaly
M = trueAnomaly2meanAnomaly_rad(ta, e);

% Return Keplerian State Vector - Mean anomaly, radians
x_kep = [sma; e; inc; raan; aop; M];

end

