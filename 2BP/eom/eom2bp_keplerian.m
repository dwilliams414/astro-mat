function dxdt = eom2bp_keplerian(t, x, GM, ad)
%EOM2BP_KEPLERIAN: Equations of motion for 2Body problem, with Keplerian
%state vector and specified disturbing acceleration
%   Inputs:
%       t:      Time of evaluation
%       x:      Keplerian state x = [a; e; i; raan; aop; M] (radians)
%       GM:     Mass parameter
%       ad:     ad(t, x, GM) Disturbing acceleration (return 3x1 vector)
   
    % Unpack state vector
   a = x(1);
   e = x(2);
   inc = x(3);
   raan = x(4);
   aop  = x(5);
   M    = x(6);

   % Compute Mean motion
   n  = sqrt(GM/x(1)^3);

   % Calculate TA and p
   ta = meanAnomaly2trueAnomaly_rad(M, e, 1e-13);
   p = a*(1-e^2);
   b = a*sqrt(1-e^2);
   h = sqrt(p*GM);

   % Compute radial magnitude
   r = p/(1+e*cos(ta));

   dxdt = zeros(6, 1);
   dxdt(end) = n;

   % Formulate B Matrix
   B = zeros(6, 3);
   
   B(1, 1) = 2*a^2*e*sin(ta);
   B(1, 2) = 2*a^2*p/r;
   B(2, 1) = p*sin(ta);
   B(2, 2) = (p+r)*cos(ta)+r*e;

   B(3, 3) = r*cos(ta+aop);
   B(4, 3) = r*sin(ta+aop)/sin(inc);

   B(5, 1) = -p*cos(ta)/e;
   B(5, 2) = (p+r)*sin(ta)/e;
   B(5, 3) = -r*sin(ta+aop)/tan(inc);

   B(6, 1) = b*p*cos(ta)/(a*e)-2*b*r/a;
   B(6, 2) = -b*(p+r)*sin(ta)/(a*e);

   B = B * 1/h;

   dxdt = dxdt + B*ad(t, x, GM);
end

