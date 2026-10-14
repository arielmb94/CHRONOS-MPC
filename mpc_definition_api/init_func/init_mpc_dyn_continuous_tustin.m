function mpc = init_mpc_dyn_continuous_tustin(mpc, Ac, Bc, Bdc, Ts)
    nx = size(Ac, 1);
    I = eye(nx);
    
    M = I - (Ts/2) * Ac;
    
    rhs_A = I + (Ts/2) * Ac;
    rhs_B = Ts * Bc;
    
    Ad = solveLinearSystemLU(M, rhs_A);
    Bd = solveLinearSystemLU(M, rhs_B);
    
    if isempty(Bdc) || (isscalar(Bdc) && Bdc == 0)
        Bdd = 0;
    else
        rhs_Bd = Ts * Bdc;
        Bdd = solveLinearSystemLU(M, rhs_Bd);
    end
    
    mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end