function phi_mats = state2phimats(states)
%STATE2PHIMATS Convert states from the result of integration to a pagewise
%array of STM's
%   Takes row-vector output from integration (states) and extracts the STM
%   at each time step, stored in a 6x6xn array
%   Inputs:
%           states: 42xn array of states from integration
%   Outputs:
%            phi_mats: 6x6xn 3D array of STM's
%       Validation: state2phimats_Validation on 09/28/2022 @ 12:48 pm
    arguments
        states (42, :) double
    end
    n_points = size(states, 2);
    stm_array = states(7:end, :);  % Need to have transpose for reshape 
                                   % to function properly
    phi_mats = reshape(stm_array, 6, 6, n_points);

end

