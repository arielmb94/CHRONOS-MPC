function mpc = update_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts)
    nx = size(Ac, 1);
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
            Ad(:,:,k) = eye(nx) + Ts * Ac(:,:,min(k, n_Ac));
            Bd(:,:,k) = Ts * Bc(:,:,min(k, n_Bc));
            if has_dist
                Bdd(:,:,k) = Ts * Bdc(:,:,min(k, n_Bdc));
            end
        end
    else
        Ad = eye(nx) + Ts * Ac;
        Bd = Ts * Bc;
        if has_dist, Bdd = Ts * Bdc; else, Bdd = 0; end
    end
    
    mpc = update_mpc_dynamics(mpc, Ad, Bd, Bdd);
end