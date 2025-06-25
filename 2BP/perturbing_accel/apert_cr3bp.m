function a3b = apert_cr3bp(t, x_pci, GM1, GM2, d_12)
%APERT_CR3BP: Model Third body perturbing acceleration, assuming a circular
%orbit of the perturbing body about the primary.
%   Inputs:
%       t:      Time of evaluation
%       x_pci:  Planet-centered inertial (PCI) spacecraft state vector
%       GM1:    Primary mass parameter
%       GM2:    Perturbing body mass parameter
%       d_12:   Circular orbit radius of perturbing body
%   Outputs:
%       a3b:        Third body perturbing acceleration, expressed in
%                   inertial (PCI) frame
    mean_motion = sqrt((GM1+GM2)/d_12^3);
    r_12 = d_12 * [cos(mean_motion*t); sin(mean_motion*t); 0];
    r = x_pci(1:3);
    
    a3b = -GM2*((r-r_12)/(norm(r-r_12)^3)+r_12/(norm(r_12)^3));
end

