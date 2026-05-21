%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [Pk_pred, x0_pred] = compute_schedul_iterative(mpc, x0, x_prev, u_prev, xref, d, sched_fun, n_rho, n_iter, A0, B0, Bd0, varargin)
%
% Computes the scheduling parameters along the prediction horizon using 
% a sequential iterative approach.
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - x0: Nx+Nu column vector, initial guess solution (warm-start).
%   - x_prev: nx column vector, current measured or estimated system state.
%   - u_prev: nu column vector, previous control action.
%   - xref: tracking reference for the MPC.
%   - d: disturbance input to the system dynamics.
%   - sched_fun: function handle (@(x)), returns n_rho scheduling parameters 
%   evaluated at state x.
%   - n_rho: integer scalar, number of scheduling variables.
%   - n_iter: integer scalar, number of sequential iterations.
%   - A0, B0, Bd0: nominal system matrices.
%   - varargin: LPV affine vertex matrices (A1..An, B1..Bn, Bd1..Bdn).
%
% Out:
%   - Pk_pred: (N * n_rho) column vector, refined scheduling trajectory.
%   - x0_pred: Nx+Nu column vector, refined state and input trajectory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pk_pred, x0_pred] = compute_schedul_iterative(mpc, x0, x_prev, u_prev, xref, d, sched_fun, ...
                                                    n_rho, n_iter, A0, B0, Bd0, varargin)
    
    % Iterative state variable (starts with provided warm-start)
    x_curr = x0;
    Pk_curr = zeros(mpc.N * n_rho, 1);
    
    % Successive trajectory refinement loop
    for iter = 1:n_iter
        
        % Extract state trajectory from current guess
        [s_curr, ~, ~] = get_x(x_curr, x_prev, mpc.nx, mpc.nu, mpc.N, mpc.N_ctr_hor, mpc.Nx);

        % TO DO: It lacks shifting
        
        % Compute Pk for current iteration based on states
        for j = 1:mpc.N
            if j == 1
                % State used to predict x(1) is the current measured state
                x_p = x_prev;
            else
                % State used to predict x(j) is the predicted x(j-1)
                idx_estado = (j-2)*mpc.nx; 
                x_p = s_curr(idx_estado + 1 : idx_estado + mpc.nx);
            end
            
            % Call user-defined scheduling function
            rho_p = sched_fun(x_p);
            
            % Populate current Pk
            Pk_curr((j-1)*n_rho + 1 : j*n_rho) = rho_p;
        end
        
        % Update dynamics with current Pk
        mpc_pred = update_mpc_sys_dynamics(mpc, A0, B0, Bd0, Pk_curr, n_rho, varargin{:});
        
        % Solve preliminary MPC to obtain refined prediction (x_next)
        [~, x_next] = mpc_solve(mpc_pred, x_curr, x_prev, u_prev, xref, d, [], [], [], ...
                                A0, Bd0, Pk_curr, n_rho, varargin{:});
                                
        % Safeguard: Abort iterations if solver fails
        % if norm(x_next) == 0
        %     Pk_pred = Pk_curr;
        %     x0_pred = x_curr;
        %     return;
        % end
        
        % Update guess for next iteration
        x_curr = x_next;
    end
    
    % Final optimized prediction after n iterations
    x0_pred = x_curr;
    
    % Extract final optimized states
    [s_final, ~, ~] = get_x(x0_pred, x_prev, mpc.nx, mpc.nu, mpc.N, mpc.N_ctr_hor, mpc.Nx);
    
    % Initialize and compute definitive scheduling vector
    Pk_pred = zeros(mpc.N * n_rho, 1);
    
    for j = 1:mpc.N
        if j == 1
            x_p = x_prev;
        else
            idx_estado = (j-2)*mpc.nx; 
            x_p = s_final(idx_estado + 1 : idx_estado + mpc.nx);
        end
        
        % Call scheduling function for final trajectory
        rho_p = sched_fun(x_p);
        
        Pk_pred((j-1)*n_rho + 1 : j*n_rho) = rho_p;
    end
end