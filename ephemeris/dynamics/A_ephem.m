function dfdx = A_ephem(true_epoch, x, model)
%A_EPHEM: Jacobian of natural ephemeris equations of motion, based upon the
%GM_nd values and instantaneous lstar.
%   Inputs:
%       true_epoch:     True epoch (time since J2000) in SECONDS
%       x:              Non-dimensional state vector
%       model:          EphemerisForceModel
%   Outputs:
%       dfdx:   Jacobian of ephemeris equations of motion
    r = x(1:3);
    rhat = r/norm(r);
    rmag = norm(r);

    ID_central = char(model.spkID(1));
    GM_central = model.GMdim(1)/model.GMstar;

    GM_nd = model.GMdim/model.GMstar;

    % Dominant Gravity Term
    dgdr = -GM_central/(rmag^3)*(eye(3)-3*(rhat*rhat'));

    % Perturbing Terms
    for k = 2:length(model.spkID)
        ID_k = char(model.spkID(k));

        r_kc = cspice_spkpos(ID_central, true_epoch, 'J2000', ...
            'None', ID_k);
        r_kc = r_kc/model.lstar; % Non-dimensionalize

        r_ks = r_kc + r;

        dgdr = dgdr - GM_nd(k)/norm(r_ks)^3*(eye(3)-3*(r_ks*r_ks')/norm(r_ks)^2);
    end

    dfdx = [zeros(3) eye(3); dgdr zeros(3)];
end