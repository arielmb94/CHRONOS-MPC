function dt_x = Furuta(x, u)
    % Nonlinear ODE of the Furuta Pendulum: Mc(x)*dt_x = f(x, u)
    [Ac, Bc, Mc] = qLPV_Furuta_SS(x);
    
    rhs = Ac * x + Bc * u;
    dt_x = solveLinearSystemLU(Mc, rhs);
end