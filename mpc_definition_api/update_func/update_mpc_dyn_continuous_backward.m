function mpc = update_mpc_dyn_continuous_backward(mpc, Ac, Bc, Bdc, Ts)
    nx = size(Ac, 1);
    I = eye(nx);
    n_Ac = size(Ac, 3);
    n_Bc = size(Bc, 3);
    N_steps = max(n_Ac, n_Bc);
    
    has_dist = ~(isempty(Bdc) || (isscalar(Bdc) && Bdc == 0));
    if has_dist
        n_Bdc = size(Bdc, 3);
    end

    if N_steps > 1
        Ad = zeros(size(Ac));
        Bd = zeros(size(Bc));
        if has_dist, Bdd = zeros(size(Bdc)); else, Bdd = 0; end
        
        for k = 1:N_steps
            M = I - Ts * Ac(:,:,min(k, n_Ac));
            
            Ad(:,:,k) = solveLinearSystemLU(M, I);
            Bd(:,:,k) = solveLinearSystemLU(M, Ts * Bc(:,:,min(k, n_Bc)));
            
            if has_dist
                Bdd(:,:,k) = solveLinearSystemLU(M, Ts * Bdc(:,:,min(k, n_Bdc)));
            end
        end
    else
        M = I - Ts * Ac;
        Ad = solveLinearSystemLU(M, I);
        Bd = solveLinearSystemLU(M, Ts * Bc);
        if has_dist, Bdd = solveLinearSystemLU(M, Ts * Bdc); else, Bdd = 0; end
    end
    
    mpc = update_mpc_dynamics(mpc, Ad, Bd, Bdd);
end