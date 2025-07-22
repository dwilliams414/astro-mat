function [M, bary_state] = rpf2inert_mat(true_epoch, model, opts)
%RPF2INERT: Determine the 6x6 matrix M that when right multipliled by a
%non-dimensional rotating-pulsating frame state, yields the corresponding
%inertial frame state, in dimensional units.  Conversion satisfies:
%   [R-B; R'-B'] = M * [rho; rho_dot]
%   Inputs:
%       true_epoch:     Epoch (seconds since J2000) at which to determine
%                       the conversion
%       model:          EphemerisForceModel.  P1 is the first ID in spkID
%                       and P2 is the second.
%   Options:
%       uniform_time:   Default is true (1).  Determines how velocity in
%                       the RPF is computed.  Non-uniform uses 
%                       t=sqrt(l^3/model.GMstar)
%   Outputs:
%       M:              6x6xn set of matrices defining conversion
%       bary_state:     6x1xn set of vectors defining barycenter state at
%                       epochs
arguments
    true_epoch (1, :) double
    model (1, 1) EphemForceModel
    opts.uniform_time (1, 1) {mustBeInteger} = 1;
end
    M = nan(6, 6, length(true_epoch));
    bary_state = nan(6, 1, length(true_epoch));
    for k = 1:size(M, 3)
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
        C = [xhat yhat zhat];

        % Angular Velocity
        %theta_dot = norm(h)/l^2;
        %omega_vec = [0; 0; theta_dot]; % angular velocity in rpf meas no.
        omega_vec = h/l^2; % express angular velocity in inertial meas no.

        % Determine M
        if opts.uniform_time
            dtdT = 1/model.tstar;
        else
            dtdT = 1/sqrt(l^3/model.GMstar);
        end
        % M(:, :, k) = [l*C zeros(3); 
        %     lprime*C+l*C*cross_matrix(omega_vec) l*C*dtdT];

        % Express DCM derivative using inertial meas no.
        M(:, :, k) = [l*C zeros(3); 
            lprime*C+l*cross_matrix(omega_vec)*C l*C*dtdT];

        % Determine Barycenter Position And Velocity
        GM2 = model.GMdim(2);

        % Barycenter state vector
        bary_state(:, :, k) = GM2/model.GMstar*x_12; 

    end
end

