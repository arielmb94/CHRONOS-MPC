function mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts, M)
    nx = size(Ac, 1);
    I = eye(nx);
    has_dist = ~(isempty(Bdc) || (isscalar(Bdc) && Bdc == 0));
    has_M = (nargin >= 6) && ~isempty(M) && ~(isscalar(M) && M == 0);
    
    if has_M
        % s+ = (I + Ts * M^(-1)*Ac)*s + Ts * M^(-1)*Bc*u + Ts * M^(-1)*Bdc*d
        Ad = I + Ts * solveLinearSystemLU(M, Ac);
        Bd = Ts * solveLinearSystemLU(M, Bc);
        if has_dist
            Bdd = Ts * solveLinearSystemLU(M, Bdc);
        else
            Bdd = 0;
        end
    else
        Ad = I + Ts * Ac;
        Bd = Ts * Bc;
        if has_dist, Bdd = Ts * Bdc; else, Bdd = 0; end
    end
    
    mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end