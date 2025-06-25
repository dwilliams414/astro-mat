function ta = eccentricAnomaly2trueAnomaly_rad(E, e)
%ECCENTRICANOMALY2TRUEANOMALY_RAD: Convert eccentric anomaly to true
%anomaly, both in radians.
%   Inputs:
%       E:      Eccentric anomaly values (radians)
%       e:      Eccentricity values (-)
%   Outputs:
%       ta:     True anomaly values (radians)
    e_term = sqrt((1+e)./(1-e));
    if (any(e >= 1))
        error("Invalid eccentricity!");
    end
    ta = 2*atan(e_term.*tan(E/2));
    ta = mod(ta, 2*pi);
end

