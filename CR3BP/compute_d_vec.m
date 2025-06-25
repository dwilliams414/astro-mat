function d = compute_d_vec(state, sys_mu)
%COMPUTE_D_VEC Vectorized code for computing CR3BP d
%   Computes the CR3BP d distance using vectorized ops, treating states as
%   row vectors
%   Validation: 09/30/2022

arguments
    state       (:, 6) double
    sys_mu      (1, 1) double
end

x = state(:, 1);
y = state(:, 2);
z = state(:, 3);

d = sqrt((x+sys_mu).^2+y.^2+z.^2);
    
end

