function aSRP = aSRP_inertial_fast(t, d_planet, GM_sun, GM_planet, AMratio, G0)
%ASRP_INERTIAL_FAST: Simplified model of solar radiation pressure, assuming
%that the Inertia n1 direction is aligned with the Sun-Planet line at time
%t = 0.
%   Inputs:
%       t:              Time since initial epoch
%       d_planet:       Radius of planet orbit about the Sun (assume
%                       a circular orbit)
%       GM_planet:      Planet mass parameter
%       AMratio:        Spacecraft area/mass ratio
%       G0:             Solar flux constant
%   Outputs:
%       aSRP:           Perturbing SRP acceleration, expressed in the 
%                       PLANET CENTERED INERTIAL FRAME!
    mean_motion = sqrt((GM_sun+GM_planet)/d_planet^3);
    s_hat = [cos(mean_motion*t); sin(mean_motion*t); 0];
    aSRP = AMratio*G0/(d_planet^2)*s_hat;
end

