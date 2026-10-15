function mpc = update_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts, M)
    nx = size(Ac, 1);
    I = eye(nx);
    
    has_dist = ~(isempty(Bdc) || (isscalar(Bdc) && Bdc == 0));
    has_M = (nargin >= 6) && ~isempty(M) && ~(isscalar(M) && M == 0);
    
    n_Ac = size(Ac, 3);
    n_Bc = size(Bc, 3);
    n_Bdc = 1; if has_dist, n_Bdc = size(Bdc, 3); end
    n_M = 1;   if has_M, n_M = size(M, 3); end
    
    N_steps = max([n_Ac, n_Bc, n_Bdc, n_M]);

    if N_steps > 1
        Ad = zeros(nx, nx, N_steps);
        Bd = zeros(nx, size(Bc, 2), N_steps);
        if has_dist, Bdd = zeros(nx, size(Bdc, 2), N_steps); else, Bdd = 0; end
        
        for k = 1:N_steps
            Ac_k = Ac(:,:,min(k, n_Ac));
            Bc_k = Bc(:,:,min(k, n_Bc));
            
            if has_M
                M_k = M(:,:,min(k, n_M));
                Ad(:,:,k) = I + Ts * solveLinearSystemLU(M_k, Ac_k);
                Bd(:,:,k) = Ts * solveLinearSystemLU(M_k, Bc_k);
                if has_dist
                    Bdc_k = Bdc(:,:,min(k, n_Bdc));
                    Bdd(:,:,k) = Ts * solveLinearSystemLU(M_k, Bdc_k);
                end
            else
                Ad(:,:,k) = I + Ts * Ac_k;
                Bd(:,:,k) = Ts * Bc_k;
                if has_dist
                    Bdd(:,:,k) = Ts * Bdc(:,:,min(k, n_Bdc));
                end
            end
        end
    else
        if has_M
            Ad = I + Ts * solveLinearSystemLU(M, Ac);
            Bd = Ts * solveLinearSystemLU(M, Bc);
            if has_dist, Bdd = Ts * solveLinearSystemLU(M, Bdc); else, Bdd = 0; end
        else
            Ad = I + Ts * Ac;
            Bd = Ts * Bc;
            if has_dist, Bdd = Ts * Bdc; else, Bdd = 0; end
        end
    end
    
    mpc = update_mpc_dynamics(mpc, Ad, Bd, Bdd);
end