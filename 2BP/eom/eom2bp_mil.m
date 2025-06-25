function dxdt = eom2bp_mil(t, x, GM, ad)
%EOM2BP_MIL: Equations of motion describing (perturbed) two-body motion
%using Milankovitch elements as the state vector, x
%   Inputs:
%       t:      Time
%       x:      Milankovitch state vector x = [h; e; L] (radians)
%       GM:     Mass parameter
%       ad:     ad(t, x, GM)->3x1 acceleration vector for perturbations
%   Outputs:
%       dxdt:   Time derivative of Milankovitch state vector
    h_vec = x(1:3);
    e_vec = x(4:6);
    L     = x(7);  % radians
    
    % Get cartesian state vector
    x_cart = mil2cart(x, GM, 1e-13); % May allow tol to change later
    r_vec = x_cart(1:3);
    v_vec = x_cart(4:6);

    % Magnitudes
    r = norm(r_vec);
    v = norm(v_vec);
    h = norm(h_vec);

    % Natural Dynamics
    f0 = zeros(7, 1);
    f0(7) = h/(r^2);

    % B - matrix
    rtilde = cross_matrix(r_vec);
    vtilde = cross_matrix(v_vec);
    htilde = cross_matrix(h_vec);


    B = [   rtilde;
            1/GM*(vtilde*rtilde-htilde);
            r_vec(3)/(h*(h+h_vec(3)))*h_vec'
        ];

    % Acceleration
    dxdt = f0 + B*ad(t, x, GM);
end