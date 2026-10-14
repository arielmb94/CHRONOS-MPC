function mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts)
    nx = size(Ac, 1);
    
    Ad = eye(nx) + Ts * Ac;
    Bd = Ts * Bc;
    
    if isempty(Bdc) || (isscalar(Bdc) && Bdc == 0)
        Bdd = 0;
    else
        Bdd = Ts * Bdc;
    end
    
    mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end