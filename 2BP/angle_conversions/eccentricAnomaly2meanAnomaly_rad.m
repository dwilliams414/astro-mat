function M = eccentricAnomaly2meanAnomaly_rad(E, e)
%ECCENTRICANOMALY2MEANANOMALY: Convert eccentric anomaly values to mean
%anomaly, both in radians. Notes G7 AAE 532
%   Inputs:
%       E:      Eccentric anomaly values (radians)
%       e:      Eccentricity values (-)
%   Outputs:
%       M:      Mean anomaly values (radians)
    if (any(e >= 1))
        error("Invalid eccentricities!");
    end
    M = E - e.*sin(E);
    M = mod(M, 2*pi);
end

