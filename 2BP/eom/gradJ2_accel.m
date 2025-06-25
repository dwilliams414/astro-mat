function dJ2 = gradJ2_accel(r, r0, GM, J2)
%GRADJ2_ACCEL Analytically construct the gradient of the J2 acceleration
%with respect to POSITION!
%   Inputs:
%       r:      Position vector
%       r0:     Primary radius
%       GM:     Primary mass parameter
%       J2:     J2 coefficient
    rm = norm(r);
    n3 = [0; 0; 1];
    I = eye(3);

    % dJ2 = rm^4*I-5*rm^2*r*r'+2*n3*n3'*(rm^4*I-5*rm^2*r*r')-...
     %    10*rm^2*(n3'*r)*r*n3'-5*(r'*n3)^2*(rm^2*I-7*r*r');

    %dJ2 = rm^4*I-5*rm^2*r*r'+2*n3*n3'*(rm^4*I-5*rm^2*r*r')-...
    %    10*rm^2*(n3'*r)*r*n3'-5*(r'*n3)^2*(rm^2*I-7*r*r');
    %dJ2 = rm^4*I-5*rm^2*r*r'+2*n3*n3'*(rm^4*I-5*rm^2*r*r')-...
    %    10*rm^2*(n3'*r)*r*n3'-5*(r'*n3)^2*(rm^2*I-7*r*r');

    dJ2 = rm^4*I-5*rm^2*r*r'+2*n3*n3'*(rm^4*I-5*rm^2*r*r')-...
        10*rm^2*(n3'*r)*r*n3'-5*(r'*n3)^2*(rm^2*I-7*r*r');
        

    dJ2 = dJ2 * -3*GM*r0^2*J2/(2*rm^9);
end

