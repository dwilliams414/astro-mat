function ad = a_valign_inert(x_inert, a_mag)
%A_VALIGN_INERT: Velocity-aligned (inertial) acceleration of fixed
%magnitude, representing a low-thrust engine
    vhat = x_inert(4:6)/norm(x_inert(4:6));
    ad = a_mag * vhat;
end