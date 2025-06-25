function M = trueAnomaly2meanAnomaly_rad(ta, e)
%TRUEANOMALY2MEANANOMALY: Convert true anomaly values to mean anomaly
%values (both radians)
%   Inputs:
%       ta:     True anomalies (radians)
%       e:      Eccentricities (-)
%   Outputs:
%       M:      Mean anomaly values
    E = trueAnomaly2eccentricAnomaly_rad(ta, e);
    M = eccentricAnomaly2meanAnomaly_rad(E, e);
end

