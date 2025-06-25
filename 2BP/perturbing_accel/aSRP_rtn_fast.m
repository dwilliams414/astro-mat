function aSRP_rtn = aSRP_rtn_fast(t, x_eci, d_planet, GM_sun, GM_planet, AMratio, G0)
%aSRP_rtn_fast: Express simplified SRP model acceleration in the RTN frame
%for use with propagating.
%   Inputs:
%       t:              Elapsed time since epoch
%       x_eci:          Spacecraft state in ECI frame (or equivalent planet-cent)
%       d_planet:       Planet orbit radius (assumed circular)
%       GM_sun:         Sun mass parameter
%       GM_planet:      Planet mass parameter
%       AMratio:        Spacecraft area to mass ratio
%   Outputs:
%       aSRP_rtn:       Solar radiation pressure acceleration in RTN frame
    aSRP_inert = aSRP_inertial_fast(t, d_planet, GM_sun, GM_planet, ...
        AMratio, G0);

    nCr = build_rtn_frame_fast(x_eci);

    aSRP_rtn = nCr'*aSRP_inert;
end