%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Pk_pred = compute_schedul_iterative(mpc, x_prev, u_prev, xref, x_ref_ter, d, sched_fun, ...
%                                       compute_A_fun, compute_B_fun, compute_Bd_fun, n_rho, n_iter)
%
% Computes the scheduling parameters along the prediction horizon using 
% a sequential iterative approach. Supports LPV matrices in A, B, and Bd.
% Includes an early stopping criterion to break iterations upon convergence.
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - x_prev: nx column vector, current measured or estimated system state.
%   - u_prev: nu column vector, previous control action.
%   - xref: tracking reference for the MPC.
%   - x_ref_ter: terminal tracking reference for the MPC.
%   - d: disturbance input to the system dynamics.
%   - sched_fun: function handle (@(x)), returns n_rho scheduling parameters 
%   evaluated at state x.
%   - compute_A_fun: function handle (@(rho)), returns matrix A. Pass [] if not LPV.
%   - compute_B_fun: function handle (@(rho)), returns matrix B. Pass [] if not LPV.
%   - compute_Bd_fun: function handle (@(rho)), returns matrix Bd. Pass [] if not LPV.
%   - n_rho: integer scalar, number of scheduling variables.
%   - n_iter: integer scalar, maximum number of sequential iterations.
%
% Out:
%   - Pk_pred: (N * n_rho) column vector, refined scheduling trajectory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pk_pred = compute_schedul_iterative(mpc, x_prev, u_prev, xref, x_ref_ter, d, sched_fun, ...
                                             compute_A_fun, compute_B_fun, compute_Bd_fun, n_rho, n_iter)
    
    % Convergence tolerance for early stopping
    tol = 1e-4;
    Pk_prev = inf(mpc.N * n_rho, 1);
    
    % Successive trajectory refinement loop
    for iter = 1:n_iter
        
        Pk_curr = zeros(mpc.N * n_rho, 1);
        
        % Compute Pk based on current mpc.s
        for j = 1:mpc.N
            if j == 1
                x_p = x_prev;
            else
                x_p = mpc.s(:, j-1);
            end
            Pk_curr((j-1)*n_rho + 1 : j*n_rho) = sched_fun(x_p);
        end
        
        % Early stopping criterion: Break if the trajectory converges
        if norm(Pk_curr - Pk_prev) < tol
            break;
        end
        Pk_prev = Pk_curr;
        
        % Build the 3D affine matrix arrays using the interface
        A_lpv = []; B_lpv = []; Bd_lpv = [];
        
        if ~isempty(compute_A_fun)
            A_lpv = traj_mat(compute_A_fun, Pk_curr, n_rho, mpc.N);
        end
        if ~isempty(compute_B_fun)
            B_lpv = traj_mat(compute_B_fun, Pk_curr, n_rho, mpc.N);
        end
        if ~isempty(compute_Bd_fun)
            Bd_lpv = traj_mat(compute_Bd_fun, Pk_curr, n_rho, mpc.N);
        end
        
        % Update dynamics with current 3D arrays
        mpc_pred = update_mpc_dynamics(mpc, A_lpv, B_lpv, Bd_lpv);
        
        % Solve preliminary MPC to obtain refined mpc structure
        % (mpc_solve now returns the updated mpc struct containing new mpc.s)
        [~, ~, mpc_next] = mpc_solve(mpc_pred, x_prev, u_prev, xref, x_ref_ter, d, [], []);
                                
        % Update mpc struct for next iteration
        mpc = mpc_next;
    end
    
    % Initialize and compute definitive scheduling vector after refinement
    Pk_pred = zeros(mpc.N * n_rho, 1);
    
    for j = 1:mpc.N
        if j == 1
            x_p = x_prev;
        else
            x_p = mpc.s(:, j-1);
        end
        Pk_pred((j-1)*n_rho + 1 : j*n_rho) = sched_fun(x_p);
    end
end