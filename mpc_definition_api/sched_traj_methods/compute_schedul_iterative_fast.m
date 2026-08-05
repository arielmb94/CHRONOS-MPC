%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pk_iterative_fast = compute_schedul_iterative_fast(mpc, x_prev, sched_fun, n_rho)
    Pk_iterative_fast = zeros(mpc.N * n_rho, 1);
    
    % Compute Pk for current iteration based on states stored in mpc.s
    for j = 1:mpc.N
        if j == 1
            % State used to predict x(1) is the current measured state
            x_p = x_prev;
        else
            % State used to predict x(j) is the predicted x(j-1)
            % Direct extraction from the state matrix mpc.s!
            x_p = mpc.s(:, j-1);
        end
        
        % Call user-defined scheduling function and populate Pk
        Pk_iterative_fast((j-1)*n_rho + 1 : j*n_rho) = sched_fun(x_p);
    end
end