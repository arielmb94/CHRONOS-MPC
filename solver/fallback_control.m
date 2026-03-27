function [u] = fallback_control(x0, x0_prev, x_prev, r, d, mpc)
    x_ref = r(1:mpc.nx);
    ref_xN = r(end-mpc.nx+1:end);
    u_incomplete = x0(1:mpc.nu);
    xN_incomplete = x0(end-mpc.nx+1:end);
    u_prev = x0_prev(mpc.nu + mpc.nx + 1: mpc.nu + mpc.nu + mpc.nx);
    xN_prev = x0_prev(end-mpc.nx+1:end);

    % Compute prediction
    pred_x_prev = mpc.A*x_prev + mpc.B*u_prev;
    pred_x_incomplete = mpc.A*x_prev + mpc.B*u_incomplete;

    if ~isempty(mpc.Bd) && ~isempty(d)
        d = d(1:mpc.nd);
        disturbance_effect = mpc.Bd*d;
        disturbance_effect = disturbance_effect(:); % for mex purposes
        pred_x_prev = pred_x_prev + disturbance_effect;
        pred_x_incomplete = pred_x_incomplete + disturbance_effect;
    end
    
    % Compute cost
    J_prev = 0;
    J_incomplete = 0;

    if ~isempty(mpc.Qe)
        J_prev = J_prev + (x_ref - pred_x_prev)'*mpc.Qe*(x_ref - pred_x_prev);
        J_incomplete = J_incomplete + (x_ref - pred_x_incomplete)'*mpc.Qe*(x_ref - pred_x_incomplete);
    end

    if ~isempty(mpc.Ru)
        J_prev = J_prev + u_prev'*mpc.Ru*u_prev;
        J_incomplete = J_incomplete + u_incomplete'*mpc.Ru*u_incomplete;
    end

    if ~isempty(mpc.ru)
        J_prev = J_prev + mpc.ru*u_prev;
        J_incomplete = J_incomplete + mpc.ru*u_incomplete;
    end

    if ~isempty(mpc.Rdu)
        J_prev = J_prev + 0;    % delta_u = 0 for u = u_prev
        J_incomplete = J_incomplete + (u_incomplete - u_prev)'*mpc.Rdu*(u_incomplete - u_prev);
    end

    if ~isempty(mpc.P)
        J_prev = J_prev + (ref_xN - xN_prev)' * mpc.P * (ref_xN - xN_prev);
        J_incomplete = J_incomplete + (ref_xN - xN_incomplete)' * mpc.P * (ref_xN - xN_incomplete);
    end

    % Evaluate smaller cost
    if J_prev < J_incomplete
        u = u_prev;
    else
        u = u_incomplete;
    end
end