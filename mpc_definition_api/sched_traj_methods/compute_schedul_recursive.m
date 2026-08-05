%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Pk_recursive = compute_schedul_recursive(mpc,x0_pred,x_curr,n_rho,sched_fun,jacob_fun,rho_min,rho_max)
%
% Computes the recursive extrapolation of qLPV scheduling parameters along 
% the prediction horizon. 
%
% Reference: 
%     Morato, M. M., Normey-Rico, J. E., & Sename, O. (2022). Sufficient 
%     conditions for convergent recursive extrapolation of qLPV scheduling 
%     parameters along a prediction horizon. IEEE Transactions on Automatic 
%     Control, 68(6), 3182-3193.
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - x0_pred: Nx+Nu column vector, predicted state and input trajectory 
%   from the previous MPC solution (used as warm-start to extract Delta x).
%   - x_curr: nx column vector, current measured or estimated system state
%   value.
%   - n_rho: integer scalar, number of scheduling variables of the qLPV 
%   model.
%   - sched_fun: function handle (@(x)), user-defined function that computes 
%   and returns the n_rho column vector of scheduling parameters rho evaluated 
%   at the state x.
%   - jacob_fun: function handle (@(x)), user-defined function that computes 
%   and returns the (n_rho x nx) Jacobian matrix (sigma_k) of the scheduling 
%   parameters evaluated at the state x.
%   - rho_min: n_rho column vector, physical lower bounds for each scheduling 
%   parameter. Used to clip the extrapolation to avoid numerical divergence.
%   - rho_max: n_rho column vector, physical upper bounds for each scheduling 
%   parameter. Used to clip the extrapolation to avoid numerical divergence.
%
% Out:
%   - Pk_recursive: (N * n_rho) column vector, extrapolated trajectory of 
%   the scheduling parameters along the prediction horizon. Formatted as 
%   [rho(k); rho(k+1); ... ; rho(k+N-1)].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pk_recursive = compute_schedul_recursive(mpc, x_curr, n_rho, sched_fun, jacob_fun, rho_min, rho_max)
    % previous optimal Pk
    prev_Pk = compute_schedul_iterative_fast(mpc, x_curr, sched_fun, n_rho);
    
    % Compute the derivative of scheduling traj. (Jacobian) sigma_k
    sigma_k = jacob_fun(x_curr); 
    
    Pk_recursive = zeros(mpc.N * n_rho, 1);
    
    % First rho is measured
    rho_curr = sched_fun(x_curr);
    Pk_recursive(1:n_rho) = rho_curr;
    
    % Recursive extrapolation along the horizon (Taylor Expansion)
    for j = 2:mpc.N
        % delta_x = x(j) - x(j-1) directly from mpc.s
        if j == 2
            x_next = mpc.s(:, 1);
            delta_x = x_next - x_curr;
        else
            x_next = mpc.s(:, j-1);
            x_prev_step = mpc.s(:, j-2);
            delta_x = x_next - x_prev_step;
        end
        
        % Extract the base rho and extrapolate
        rho_base = prev_Pk((j-1)*n_rho + 1 : j*n_rho);
        rho_next = rho_base + sigma_k * delta_x;
        
        % Clip each entry of the extrapolation
        rho_next = max(rho_min, min(rho_next, rho_max));
        
        % Complete Pk
        idx_rho = (j-1)*n_rho;
        Pk_recursive(idx_rho + 1 : idx_rho + n_rho) = rho_next;
    end
end