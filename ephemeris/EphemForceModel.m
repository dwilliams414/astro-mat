classdef EphemForceModel < handle
    %EphemForceModel:  Class for describing parameters associated with
    %ephemeris EOM's and propagation.  This is predominantly geared toward
    %uniform time flow in the inertial frame, and modifications will be
    %necessary when the RPS implementation is incorporated.

    properties
        spkID (1, :) string;        % SPK ID values
        GMdim (1, :) double;        % GM values for bodies (dimensional)
        epoch_str (1, 1) string;    % Initial epoch as datetime string
        et0   (1, 1) double;        % Ephemeris time for initial epoch
        tstar (1, 1) double;        % Characteristic time (fixed val)
        lstar (1, 1) double;        % Characteristic length (fixed val)
        GMstar (1, 1) double;       % Characteristic GM value (fixed val)
    end

    methods
        function obj = EphemForceModel(opts)
            %EphemForceModel Construct an instance of this class
            %   Options:
            %       'spkID':        Array of strings giving the SPK ID
            %                       values of bodies used in the model
            %
            %       'GMdim':        Array of double describing the mass
            %                       parameters of bodies used in model.
            %                       This overwrites the default SPICE
            %                       values, which is generally a bad
            %                       practice.  Given in [km^3/s^2]
            %       
            %       'epoch_str':    String defining initial epoch
            %
            %       'et0':          Epoch time (seconds since J2000)
            %
            %       
            arguments
                opts.spkID (1, :) string = ["Earth", "Moon", "Sun"];
                opts.GMdim (1, :) double = [];
                opts.epoch_str (1, 1) string = "20 Apr 2025 00:00:00.000";
                opts.et0 (1, 1) double = inf;
                opts.lstar (1, 1) double = inf;
            end

            % Configure bodies for use in model
            obj.spkID = opts.spkID;

            % Get Dimensional GM Values
            if isempty(opts.GMdim)
                % Set default values from SPICE data
                obj.GMdim = zeros(size(opts.GMdim));
                for k = 1:length(obj.spkID)
                    obj.GMdim(k) = cspice_bodvrd(char(obj.spkID(k)), 'GM', 1);
                end
            else
                obj.GMdim = opts.GMdim;
            end

            % Calculate Mass Ratio/GMstar
            if length(obj.spkID) > 1
                obj.GMstar = (obj.GMdim(1)+obj.GMdim(2));
            else
                obj.GMstar = obj.GMdim(1);
            end

            % Determine Initial Epoch
            if isinf(opts.et0)
                obj.epoch_str = opts.epoch_str;
                obj.et0 = cspice_str2et(char(obj.epoch_str));
            else
                obj.et0 = opts.et0;
                obj.epoch_str = string(cspice_et2utc(obj.et0, 'C', 4));
            end

            % Determine lstar (for propagation, fixed val for now)
            if isinf(opts.lstar)
                r_12 = cspice_spkpos(char(obj.spkID(2)), obj.et0, ...
                    'J2000', 'None', char(obj.spkID(1)));
                obj.lstar = norm(r_12);
            else
                obj.lstar = opts.lstar;
            end

            % Determine tstar
            obj.tstar = sqrt(obj.lstar^3/obj.GMstar);
        end

        function x_dim = dimensionalize(obj, x_nd)
            arguments
                obj  (1, 1) EphemForceModel
                x_nd (6, :)
            end
            dim_mat = [obj.lstar*eye(3) zeros(3);
                zeros(3) obj.lstar*eye(3)/obj.tstar];

            x_dim = dim_mat*x_nd;
        end

        function x_nd = nondimensionalize(obj, x_dim)
            arguments
                obj(1, 1) EphemForceModel
                x_dim (6, :) double
            end
            nd_mat = [1/obj.lstar*eye(3) zeros(3);
                zeros(3) obj.tstar/obj.lstar*eye(3)];
            x_nd = nd_mat*x_dim;
        end

        function t_true = get_true_epochs(obj, t_nd)
           arguments
                obj (1, 1) EphemForceModel
                t_nd (1, :) double
           end
           t_true = obj.et0 + t_nd*obj.tstar;
        end

        function t_nd = get_nd_time_since_epoch(obj, true_epochs)
            arguments
                obj (1, 1) EphemForceModel
                true_epochs (1, :) double
            end
            t_nd = (true_epochs-obj.et0)/obj.tstar;
        end
    end
end