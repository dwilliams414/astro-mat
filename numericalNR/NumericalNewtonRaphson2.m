classdef NumericalNewtonRaphson2 < handle
    %NUMERICALNEWTONRAPHSON Class for fully numerical solution of
    %Newton-Raphson problem with NO analytical gradients.  Constructor
    %requires only a function handle, F(X), where X is some vector.
    
    properties
        F_fun; % F(X) function to solve (X is a vector)
        h;     % Step Size (or vector of step-sizes)
        max_iter; % Maximum number of iterations for solving
        err_tol;  % Norm error tolerance OR vector of individual tolerances
        use_comp_error; % Flag to use or not use component errors
        error_on_fail; % Throw an error if solution does not converge
        use_forward_diff; % Use forward as opposed to central difference
        use_max_update;   % Use a maximum update size (norm)
        max_update_size;  % Maximum allowable size for update
        update_factor;    % Scale factor to multiply by if update too large
                          % and using relative update (as opposed to abs)
        use_rel_update;   % Flag to use relative as opposed to absolute 
                          % update resizing in solve().
    end
    
    methods
        % Constructor
        function obj = NumericalNewtonRaphson2(F_fun, options)
            arguments
                F_fun function_handle;
                options.max_iter = 25;
                options.h        = sqrt(eps);
                options.err_tol (:, 1)  = 1e-12;
                options.error_on_fail = 1;
                options.use_forward_diff = 0;
                options.max_update_size = inf;
                options.use_max_update  = 0;
                options.update_factor = 0.5;
                options.use_rel_update (1, 1) {mustBeInteger} = 0;
            end
            %NUMERICALNEWTONRAPHSON Construct an instance of this class
            obj.F_fun = F_fun;
            obj.h = options.h;
            obj.max_iter = options.max_iter;
            obj.err_tol = options.err_tol;
            obj.error_on_fail = options.error_on_fail;
            obj.use_max_update = options.use_max_update;
            obj.max_update_size = options.max_update_size;
            obj.update_factor   = options.update_factor;
            obj.use_rel_update = options.use_rel_update;

            
            if (length(options.err_tol) > 1)
                obj.use_comp_error = 1;
            else
                obj.use_comp_error = 0;
            end
        end
        
        % Finite Differencing Utility
        function DF = finite_diff(obj,X0)

            %finite_diff: Generate finite differenced DF matrix for X0
            arguments
                obj NumericalNewtonRaphson2
                X0 (:, 1) double
            end

            % Size Matrices
            F_nominal = obj.F_fun(X0);
            DF = zeros(size(F_nominal, 1), length(X0));

            % Apply Finite Differencing
            for ii = 1:length(X0)
                X0p = X0;
                X0m = X0;

                if (length(obj.h) == 1)
                    X0p(ii) = X0(ii)+obj.h;
                    X0m(ii) = X0(ii)-obj.h;
                elseif (length(obj.h) == length(X0))
                    X0p(ii) = X0(ii)+obj.h(ii);
                    X0m(ii) = X0(ii)-obj.h(ii);
                else
                    error("Invalid h specified in NumericalNewtonRaphson.");
                end

                Fp = obj.F_fun(X0p);
                Fm = obj.F_fun(X0m);

                if (length(obj.h) == 1)
                    DF(:, ii) = (Fp-Fm)/(2*obj.h);
                else
                    DF(:, ii) = (Fp-Fm)/(2*obj.h(ii));
                end
            end


        end
        
        % Solver
        function [X_solved, F_solved, iter_error] = solve(obj, X0)
            arguments
                obj NumericalNewtonRaphson2
                X0 (:, 1) double
            end

            % Variables
            X_i = X0; % Current Iteration Guess

            % Evaluate F for Initial Guess
            F_i = obj.F_fun(X0);

            % Ensure appropriate error vector has been specified
            if (obj.use_comp_error && length(F_i) ~= length(obj.err_tol))
                warning("Dimension mismatch for component-wise error" + ...
                    "checking.  Switching to Norm");
                obj.use_comp_error = 0;
            end

            % Initialize num_iter and iter_error
            num_iter = 0;

            if (~obj.use_comp_error)
                iter_error = nan(obj.max_iter, 1);
                iter_error(num_iter+1) = norm(F_i);
            else
                iter_error = nan(obj.max_iter, length(F_i));
                iter_error(num_iter+1, :) = F_i;
            end

            % Iterate
            while ((norm(F_i) > obj.err_tol(1) && ~obj.use_comp_error) ...
                    || (any(abs(F_i) > obj.err_tol) && obj.use_comp_error) && ...
                    (num_iter < obj.max_iter))

                if obj.use_forward_diff
                    DF = obj.forward_diff(X_i);
                else
                    DF = obj.finite_diff(X_i);
                end

                if (size(DF, 1) == size(DF, 2))
                    update = -DF\F_i;
                else
                    update = lsqminnorm(-DF, F_i);
                end

                % Check if update is to be reduced
                while (obj.use_max_update && norm(update) > obj.max_update_size+eps())
                    if (obj.use_rel_update)
                        update = update * obj.update_factor;
                    else
                        update = update/norm(update) * obj.max_update_size;
                    end
                end

                % Update Estimate and assoc. F
                X_i = X_i + update;
                F_i = obj.F_fun(X_i);

                % Increment Iterations & iter_error
                num_iter = num_iter + 1;

                if (~obj.use_comp_error)
                    iter_error(num_iter+1) = norm(F_i);
                else
                    iter_error(num_iter+1, :) = F_i;
                end
            end

            if (num_iter < obj.max_iter)
                X_solved = X_i;
                F_solved = F_i;
            elseif (obj.error_on_fail)
                error("Solution Failed to Converge!");
            else
                X_solved = nan(size(X0));
                F_solved = nan(size(F_i));
                warning("Solution Failed to Converge!  Returning nan");
            end
        end

        % Forward Differencing Utility
        function DF = forward_diff(obj, X0)
            arguments
                obj NumericalNewtonRaphson2
                X0 (:, 1) double
            end

            % Size Matrices
            F_nominal = obj.F_fun(X0);
            DF = zeros(size(F_nominal, 1), length(X0));

            % Apply Finite Differencing
            for ii = 1:length(X0)
                X0p = X0;

                if (length(obj.h) == 1)
                    X0p(ii) = X0(ii)+obj.h;
                elseif (length(obj.h) == length(X0))
                    X0p(ii) = X0(ii)+obj.h(ii);
                else
                    error("Invalid h specified in NumericalNewtonRaphson.");
                end

                Fp = obj.F_fun(X0p);

                if (length(obj.h) == 1)
                    DF(:, ii) = (Fp-F_nominal)/(obj.h);
                else
                    DF(:, ii) = (Fp-F_nominal)/(obj.h(ii));
                end
            end % for
        end % function
    end % methods
end % classdef

