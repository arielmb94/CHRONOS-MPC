function mpc = f0_hess_tracking_wrapper(mpc)
[mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Q_f0_ter] = f0_hess_tracking_local(mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Q_f0_ter,mpc.y_use_k0,mpc.y_use_s,mpc.y_use_u,mpc.y_use_ter,mpc.R_E_0,mpc.Q_E,mpc.R_E,mpc.Y_E,mpc.Q_E_ter,mpc.s_col);
end

function [R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter] = f0_hess_tracking_local(R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter,y_use_k0,y_use_s,y_use_u,y_use_ter,R_E_0,Q_E,R_E,Y_E,Q_E_ter,s_col)
if ~isempty(y_use_k0)
    R_f0_0(:,:) = R_f0_0 + R_E_0;
end
if ~isempty(y_use_s) && ~isempty(y_use_u)
    Q_f0_k(s_col,s_col,:) = Q_f0_k(s_col,s_col,:) + Q_E;
    R_f0_k(:,:,:) = R_f0_k + R_E;
    Y_f0_k(:,s_col,:) = Y_f0_k(:,s_col,:) + Y_E;
elseif ~isempty(y_use_s)
    Q_f0_k(s_col,s_col,:) = Q_f0_k(s_col,s_col,:) + Q_E;
elseif ~isempty(y_use_u)
    R_f0_k(:,:,:) = R_f0_k + R_E;
end
if ~isempty(y_use_ter)
    Q_f0_ter(s_col,s_col) = Q_f0_ter(s_col,s_col) + Q_E_ter;
end
end
