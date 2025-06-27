function dxdt = f_ephem_uni_natural(t, x, model)
%f_ephem: Natural ephemeris equations of motion, formulated in the inertial
%J2000 frame and assuming UNIFORM TIME FLOW!
%Note that this formulation has no input argument validation, and generally
%should not be called by the user UNLESS they know what they're doing.  All
%propagation occurs in the J2000 frame.
%   Inputs:
%       t:          Non-dimensionalized time since epoch (et)
%       x:          Non-dimensionalized state vector (6x1, 42x1, or 48x1)
%       model:      EphemForceModel
%   Outputs:
%       dxdt:       State derivative

    r = x(1:3);
    rmag = norm(r);

    % Unpack Variables
    epoch_time = model.et0;
    GM_dim = model.GMdim;
    spk_ID = model.spkID;
    tstar = model.tstar;
    lstar = model.lstar;

    % Allocate for Derivative
    dxdt = zeros(length(x), 1);
    dxdt(1:3) = x(4:6);
    % dxdt = dxdt';

    % Dominant Term
    GM_ndim = GM_dim / lstar^3 * tstar^2;
    dxdt(4:6) = -GM_ndim(1)/rmag^3*r;
    ID_c = char(spk_ID(1));

    % Perturbing
    for k = 2:length(spk_ID)
        ID_k = char(spk_ID(k));
        r_kc = cspice_spkpos(ID_c, epoch_time+t*tstar, 'J2000', 'None', ...
            ID_k);

        r_kc = r_kc/lstar;
        r_ks = r_kc+r;

        dxdt(4:6) = dxdt(4:6)-GM_ndim(k)*(r_ks/(norm(r_ks)^3)-r_kc/(norm(r_kc)^3));
    end

    % STM Elements (if required)
    if length(x) > 6
        A = A_ephem(epoch_time+tstar*t, x(1:6), model);
        phi = reshape(x(7:42), 6, 6);
        phidot = A*phi;
        dxdt(7:42) = phidot(:);
    end
    
end