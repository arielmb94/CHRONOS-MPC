function [u] = fallback_control(u_prev,u_unconstrained, x_prev, x_ref, mpc)
    pred_x_1 = mpc.A*x_prev + mpc.B*u_prev;
    pred_x_2 = mpc.A*x_prev + mpc.B*u_unconstrained(1);
    J_1 = (x_ref - pred_x_1)'*mpc.Qe*(x_ref - pred_x_1);
    J_2 = (x_ref - pred_x_2)'*mpc.Qe*(x_ref - pred_x_2);

    if J_1 < J_2
        u = u_prev;
    else
        u = u_unconstrained;
    end
end