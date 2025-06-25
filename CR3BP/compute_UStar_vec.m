function UStar = compute_UStar_vec(state, sys_mu)
%COMPUTE_USTAR_VEC Vectorized computation of CR3BP pseudopotential

arguments
    state (:, 6) double
    sys_mu (1, 1) double
end
d = compute_d_vec(state, sys_mu);
r = compute_r_vec(state, sys_mu);

UStar = (1-sys_mu)./d + sys_mu./r + 0.5*(state(:,1).^2+state(:,2).^2);

end

