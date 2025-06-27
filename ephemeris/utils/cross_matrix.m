function v_tilde = cross_matrix(vec)
%CROSS_MATRIX: Return the matrix representation of a cross product.  I.E.,
%if vec is a vector, return rtilde such that rvec x b = rtilde*b
arguments
    vec(3, :) double
end
    v_tilde = zeros(3, 3, size(vec, 2));

    v_tilde(1, 2, :) = -vec(3, :);
    v_tilde(1, 3, :) = vec(2, :);

    v_tilde(2, 1, :) = vec(3, :);
    v_tilde(2, 3, :) = -vec(1, :);

    v_tilde(3, 1, :) = -vec(2, :);
    v_tilde(3, 2, :) = vec(1, :);
end

