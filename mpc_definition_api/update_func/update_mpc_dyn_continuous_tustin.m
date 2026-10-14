function mpc = update_mpc_dyn_continuous_tustin(mpc, Ac, Bc, Bdc, Ts)
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
            Ac_k = Ac(:,:,min(k, n_Ac));
            M = I - (Ts/2) * Ac_k;
            
            rhs_A = I + (Ts/2) * Ac_k;
            rhs_B = Ts * Bc(:,:,min(k, n_Bc));
            
            Ad(:,:,k) = solveLinearSystemLU(M, rhs_A);
            Bd(:,:,k) = solveLinearSystemLU(M, rhs_B);
            
            if has_dist
                rhs_Bd = Ts * Bdc(:,:,min(k, n_Bdc));
                Bdd(:,:,k) = solveLinearSystemLU(M, rhs_Bd);
            end
        end
    else
        M = I - (Ts/2) * Ac;
        rhs_A = I + (Ts/2) * Ac;
        
        Ad = solveLinearSystemLU(M, rhs_A);
        Bd = solveLinearSystemLU(M, Ts * Bc);
        if has_dist, Bdd = solveLinearSystemLU(M, Ts * Bdc); else, Bdd = 0; end
    end
    
    mpc = update_mpc_dynamics(mpc, Ad, Bd, Bdd);
end