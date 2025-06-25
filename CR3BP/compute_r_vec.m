function r = compute_r_vec(states,sys_mu)
%COMPUTE_R_VEC Vectorized function for computing "r" distance in CR3BP
%       Validation: 09/30/2022
arguments
    states (:, 6) double
    sys_mu (1, 1) double
end

x = states(:, 1);
y = states(:, 2);
z = states(:, 3);

r = sqrt((x+sys_mu-1).^2+y.^2+z.^2);

end