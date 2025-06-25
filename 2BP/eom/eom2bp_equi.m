function dxdt = eom2bp_equi(t, x, GM, ad)
%EOM2BP_EQUI: Two body problem equations of motion in modified equinoctial
%orbital elements, with a perturbing acceleration
%   Inputs:
%       t:      Time of evaluation
%       x:      Modified Equinoctial state vector
%               x = [p f g h k L]' (L in radians!)
%       GM:     Mass paramter of primary
%       ad:     ad(t, x, GM)-> 3x1 vector defining the disturbing
%               acceleration
%   Outputs:
%       dxdt:   Modified equinoctial state vector time derivative

    % Unpack state elements
    p = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    L = x(6);

    % Addl Parameters
    q = 1 + f*cos(L) + g*sin(L);
    ssq = 1 + h^2 + k^2;

    % Natural Dynamics
    f0 = zeros(6, 1);
    f0(6) = sqrt(GM*p)*(q/p)^2;

    % B Matrix
    B = zeros(6, 3);
    B(1, 2) = 2*p/q;

    B(2, 1) = sin(L);
    B(2, 2) = ((q+1)*cos(L)+f)/q;
    B(2, 3) = -(g*(h*sin(L)-k*cos(L)))/q;

    B(3, 1) = -cos(L);
    B(3, 2) = ((q + 1)*sin(L) + g)/q;
    B(3, 3) = f*(h*sin(L) - k*cos(L))/q;

    B(4, 3) = ssq/(2*q)*cos(L);
    B(5, 3) = ssq/(2*q)*sin(L);
    B(6, 3) = (h*sin(L)-k*cos(L))/q;

    dxdt = f0 + sqrt(p/GM)*B*ad(t, x, GM);

end

