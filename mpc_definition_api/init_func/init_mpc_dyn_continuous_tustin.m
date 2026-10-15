function mpc = init_mpc_dyn_continuous_tustin(mpc, Ac, Bc, Bdc, Ts, M)
    nx = size(Ac, 1);
    has_dist = ~(isempty(Bdc) || (isscalar(Bdc) && Bdc == 0));
    if (nargin < 6) || isempty(M) || (isscalar(M) && M == 0)
        M_mat = eye(nx);
    else
        M_mat = M;
    end
    
    % (M - (Ts/2)*Ac)*s+ = (M + (Ts/2)*Ac)*s + Ts*Bc*u + Ts*Bdc*d
    W = M_mat - (Ts/2) * Ac;
    rhs_A = M_mat + (Ts/2) * Ac;
    
    Ad = solveLinearSystemLU(W, rhs_A);
    Bd = solveLinearSystemLU(W, Ts * Bc);
    if has_dist
        Bdd = solveLinearSystemLU(W, Ts * Bdc);
    else
        Bdd = 0;
    end
    
    mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end