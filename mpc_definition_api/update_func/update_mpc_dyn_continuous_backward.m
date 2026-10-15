function mpc = update_mpc_dyn_continuous_backward(mpc, Ac, Bc, Bdc, Ts, M)
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
            if has_M, M_k = M(:,:,min(k, n_M)); else, M_k = I; end
            
            W = M_k - Ts * Ac_k;
            Ad(:,:,k) = solveLinearSystemLU(W, M_k);
            Bd(:,:,k) = solveLinearSystemLU(W, Ts * Bc_k);
            if has_dist
                Bdc_k = Bdc(:,:,min(k, n_Bdc));
                Bdd(:,:,k) = solveLinearSystemLU(W, Ts * Bdc_k);
            end
        end
    else
        if has_M, M_mat = M; else, M_mat = I; end
        W = M_mat - Ts * Ac;
        
        Ad = solveLinearSystemLU(W, M_mat);
        Bd = solveLinearSystemLU(W, Ts * Bc);
        if has_dist, Bdd = solveLinearSystemLU(W, Ts * Bdc); else, Bdd = 0; end
    end
    
    mpc = update_mpc_dynamics(mpc, Ad, Bd, Bdd);
end