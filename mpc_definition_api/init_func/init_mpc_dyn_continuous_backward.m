function mpc = init_mpc_dyn_continuous_backward(mpc, Ac, Bc, Bdc, Ts)
    nx = size(Ac, 1);
    I = eye(nx);
    
    M = I - Ts * Ac;
    
    % Solve MX = B using LU decomposition
    Ad = solveLinearSystemLU(M, I);
    Bd = solveLinearSystemLU(M, Ts * Bc);
    
    if isempty(Bdc) || (isscalar(Bdc) && Bdc == 0)
        Bdd = 0;
    else
        Bdd = solveLinearSystemLU(M, Ts * Bdc);
    end
    
    mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end