classdef CR3BPSystem
    %CR3BPSystem: CR3BP System definition
    %   Describes a CR3BP system.

    properties
        lstar;
        tstar;
        mass_ratio;
        GMtotal;
    end

    methods
        function obj = CR3BPSystem(lstar, GMtotal, mass_ratio)
            % Inputs:
            %       lstar: Characteristic length
            %       GMtotal: Total mass parameter G(M1+M2) for primaries
            %       mass_ratio: GM2/(GM1+GM2)
            % Outputs:
            %       Instance of CR3BP system
            arguments
                lstar (1, 1) double
                GMtotal (1, 1) double
                mass_ratio (1, 1) double
            end
            obj.lstar = lstar;
            obj.GMtotal = GMtotal;
            obj.lstar = lstar;
            obj.mass_ratio = mass_ratio;
            obj.tstar = 1/sqrt(GMtotal/lstar^3);

        end

        function x_dim = dimensionalize_states(obj, x_nd)
            %Dimensionalize states according to the characteristic length
            %and time associated with this system
            %   Inputs:
            %       x_nd:   Nondimensional states
            %   Outputs:
            %       x_dim:  Dimensionalized States
            arguments
                obj  (1, 1) CR3BPSystem
                x_nd (6, :) double
            end
            x_dim = x_nd;
            x_dim(1:3, :) = x_dim(1:3, :)*obj.lstar;

            x_dim(4:6, :) = x_dim(4:6, :)*obj.lstar/obj.tstar;

        end

        function x_nd = nondimensionalize_states(obj, x_dim)
            %Dimensionalize states according to the characteristic length
            %and time associated with this system
            %   Inputs:
            %       x_dim:      Dimensional states
            %   Outputs:
            %       x_nd:      Nondim States
            arguments
                obj  (1, 1) CR3BPSystem
                x_dim (6, :) double
            end
            
            x_nd = x_dim;
            x_nd(1:3, :) = x_dim(1:3, :)/obj.lstar;
            x_nd(4:6, :) = x_dim(4:6, :)/obj.lstar*obj.tstar;

        end
    end
end