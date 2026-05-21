%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pk_iterative_fast = compute_schedul_iterative_fast(mpc, x0, x_prev, sched_fun, n_rho)
    x_curr = x0;
    Pk_iterative_fast = zeros(mpc.N * n_rho, 1);
    
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
        
        % Populate Pk
        Pk_iterative_fast((j-1)*n_rho + 1 : j*n_rho) = rho_p;
    end
end