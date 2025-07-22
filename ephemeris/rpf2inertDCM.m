function nCr = rpf2inertDCM(true_epoch, model)
%rpf2inertDCM: Build the instantaneous DCM between the J2000 and RPF
%   Inputs:
%       true_epoch:     Epoch of conversion (seconds since J2000)
%       model:          EphemForceModel
%   Outputs:
%       nCr:            [xhat, yhat, zhat] defining the RPF
   arguments
        true_epoch (1, 1) double
        model (1, 1) EphemForceModel
   end
    % Determine characteristic length and derivative wrt dimensional
    % time
    x_12 = cspice_spkezr(char(model.spkID(2)), true_epoch(k), ...
        'J2000', 'None', char(model.spkID(1)));
    r12 = x_12(1:3);
    v12 = x_12(4:end);

    % Characteristic Length and Derivative
    l = norm(r12);
    lprime = (r12/l)'*v12;
    xhat = r12/norm(r12);

    % Determine Angular Momentum Vector (unit)
    h = cross(r12, v12);
    zhat = h/norm(h);

    % Complete Triad
    yhat = cross(zhat, xhat);

    % DCM
    nCr = [xhat yhat zhat];
end