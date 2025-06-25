function ta = meanAnomaly2trueAnomaly_rad(M, e, tol)
%MEANANOMALY2TRUEANOMALY_RAD: Convert Mean anomaly values to true anomaly
%values (both radians):
%   Inputs:
%       M:          Mean anomaly values (rad)
%       e:          Eccentricity values (-)
%   Outputs:
%       ta:         True anomaly values (radians)
    E = meanAnomaly2eccentricAnomaly_rad(M, e, tol);
    ta = eccentricAnomaly2trueAnomaly_rad(E, e);
end

