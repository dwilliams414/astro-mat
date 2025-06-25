function dcm_313 = body2_313(theta1, theta2, theta3)
%DCM_313: Get the direction cosine matrix for the 313 rotation sequence.
%Angles are given in RADIANS!
%   Inputs:
%       theta1:     First rotation angle (rad)
%       theta2:     Second rotation angle (rad)
%       theta3:     Third rotation angle (rad)
% arguments
%     theta1 (:, 1) double
%     theta2 (:, 1) double
%     theta3 (:, 1) double {validateDimension(theta1, theta2, theta3)}
% end
    dcm_313 = zeros(3, 3, length(theta1));
    t1 = reshape(theta1, 1, 1, []);
    t2 = reshape(theta2, 1, 1, []);
    t3 = reshape(theta3, 1, 1, []);

    % First Row
    dcm_313(1, 1, :) = -sin(t1).*cos(t2).*sin(t3)+cos(t3).*cos(t1);
    dcm_313(1, 2, :) = -sin(t1).*cos(t2).*cos(t3)-sin(t3).*cos(t1);
    dcm_313(1, 3, :) = sin(t1).*sin(t2);

    dcm_313(2, 1, :) = cos(t1).*cos(t2).*sin(t3)+cos(t3).*sin(t1);
    dcm_313(2, 2, :) = cos(t1).*cos(t2).*cos(t3)-sin(t3).*sin(t1);
    dcm_313(2, 3, :) = -cos(t1).*sin(t2);

    dcm_313(3, 1, :) = sin(t2).*sin(t3);
    dcm_313(3, 2, :) = sin(t2).*cos(t3);
    dcm_313(3, 3, :) = cos(t2);

%     % Initial Implementation
%     rotation_1 = [cos(theta1) -sin(theta1) 0;
%                   sin(theta1) cos(theta1) 0;
%                   0 0 1];
% 
%     rotation_2 = [1 0 0;
%                   0 cos(theta2) -sin(theta2);
%                   0 sin(theta2) cos(theta2)];
%     rotation_3 = [cos(theta3) -sin(theta3) 0;
%                   sin(theta3) cos(theta3) 0;
%                   0 0 1];
% 
%     dcm_313 = rotation_1 * rotation_2 * rotation_3;
end

function validateDimension(theta1, theta2, theta3)
    num_dcm = length(theta1);

    if (length(theta2)~=num_dcm || length(theta3)~=num_dcm || ...
            length(theta2)~=length(theta3))
        error("Invalid angles specified!  Dimensions much match");
    end
end

