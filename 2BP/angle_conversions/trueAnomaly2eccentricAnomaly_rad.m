function E = trueAnomaly2eccentricAnomaly_rad(ta, ecc)
%TRUEANOMALY2ECCENTRICANOMAL: Convert true anomaly to eccentric anomaly
%(both in radians).  Based on formula from notes G2 AAE 532
%   Inputs:
%       ta:     True anomaly values (radians)
%       ecc:    Eccentricity values (-)
%   Outputs
%       E:      Eccentric anomaly values (radians)
    if (any(ecc >= 1))
        error("Invalid eccentricities!");
    end
    ecc_term = sqrt((1+ecc)./(1-ecc));
    E = 2*atan(1./ecc_term.*tan(ta/2));
    E = mod(E, 2*pi);
end

