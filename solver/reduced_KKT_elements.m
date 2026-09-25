function mpc = reduced_KKT_elements(mpc,has_u_cnstr,has_du_cnstr,has_s_cnstr,has_y_cnstr,has_h_cnstr,g_0,g_k,g_ter)

[mpc.ru_0,mpc.rse_k,mpc.ru_k,mpc.rse_ter,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter] = ...
    reduced_kkt_objective_local(mpc.t,mpc.R_0,mpc.Q_k,mpc.R_k,mpc.Y_k,mpc.Q_ter,...
    mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Q_f0_ter,mpc.ru_0,...
    mpc.rse_k,mpc.ru_k,mpc.rse_ter,mpc.grad_u_f0_0,mpc.grad_se_f0_k,...
    mpc.grad_u_f0_k,mpc.grad_se_f0_ter);

if ~isempty(g_0)
    [mpc.g2_0,mpc.v2_0,mpc.rv_v2_0,mpc.ri_0,mpc.iS_0,mpc.iS_ri_hat_0] = ...
        reduced_kkt_0_local(mpc.t,mpc.g_0,mpc.v_0,mpc.g2_0,mpc.v2_0,...
        mpc.rv_v2_0,mpc.ri_0,mpc.grad_qv_0,mpc.v_rows_0,mpc.iS_0,...
        mpc.iS_ri_hat_0);
end
if ~isempty(g_k)
    [mpc.g2_k,mpc.v2_k,mpc.rv_v2_k,mpc.ri_k,mpc.iS_k,mpc.iS_ri_hat_k] = ...
        reduced_kkt_k_local(mpc.t,mpc.g_k,mpc.v_k,mpc.g2_k,mpc.v2_k,...
        mpc.rv_v2_k,mpc.ri_k,mpc.grad_qv_k,mpc.v_rows_k,mpc.iS_k,...
        mpc.iS_ri_hat_k);
end
if ~isempty(g_ter)
    [mpc.g2_ter,mpc.v2_ter,mpc.rv_v2_ter,mpc.ri_ter,mpc.iS_ter,...
        mpc.iS_ri_hat_ter] = reduced_kkt_terminal_local(mpc.t,mpc.g_ter,...
        mpc.v_ter,mpc.g2_ter,mpc.v2_ter,mpc.rv_v2_ter,mpc.ri_ter,...
        mpc.grad_qv_ter,mpc.iS_ter,mpc.iS_ri_hat_ter);
end
if ~isempty(has_u_cnstr), mpc = reduced_kkt_control_wrapper(mpc,mpc.u_cnstr); end
if ~isempty(has_du_cnstr), mpc = reduced_kkt_control_rate_wrapper(mpc,mpc.du_cnstr); end
if ~isempty(has_s_cnstr), mpc = reduced_kkt_state_wrapper(mpc,mpc.s_cnstr); end
if ~isempty(has_y_cnstr), mpc = reduced_kkt_output_wrapper(mpc,mpc.y_cnstr); end
if ~isempty(has_h_cnstr), mpc = reduced_kkt_custom_wrapper(mpc,mpc.h_cnstr); end
end

function [ru_0,rse_k,ru_k,rse_ter,R_0,Q_k,R_k,Y_k,Q_ter] = ...
    reduced_kkt_objective_local(t,R_0,Q_k,R_k,Y_k,Q_ter,R_f0_0,Q_f0_k,...
    R_f0_k,Y_f0_k,Q_f0_ter,ru_0,rse_k,ru_k,rse_ter,grad_u_f0_0,...
    grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter)
ru_0 = grad_u_f0_0;
rse_k = grad_se_f0_k;
ru_k = grad_u_f0_k;
rse_ter = grad_se_f0_ter;

R_0 = t * R_f0_0;
Q_k = t * Q_f0_k;
R_k = t * R_f0_k;
Y_k = t * Y_f0_k;
Q_ter = t * Q_f0_ter;
end

function [g2_0,v2_0,rv_v2_0,ri_0,iS_0,iS_ri_hat_0] = ...
    reduced_kkt_0_local(t,g_0,v_0,g2_0,v2_0,rv_v2_0,ri_0,grad_qv_0,...
    v_rows_0,iS_0,iS_ri_hat_0)
%% ri_hat = ri - g^2*(-1/g) + v^2*(t*qv-1/v)
%  ri_hat = ri + g + t*qv*v^2 - v

ri_0 = ri_0+g_0;
g2_0 = g_0.^2;
if ~isempty(v_0)
    v2_0 = v_0.^2;
    % v^2*rv = t*qv*v^2 - v
    rv_v2_0 = t*grad_qv_0.*v2_0-v_0;
    ri_0(v_rows_0) = ri_0(v_rows_0) + rv_v2_0;
end

%% S = g^2 + v^2
iS_0 = g2_0;
if ~isempty(v_0), iS_0(v_rows_0) = iS_0(v_rows_0) + v2_0; end
iS_0  = 1./iS_0;
iS_ri_hat_0 = iS_0.*ri_0;
end

function [g2_k,v2_k,rv_v2_k,ri_k,iS_k,iS_ri_hat_k] = ...
    reduced_kkt_k_local(t,g_k,v_k,g2_k,v2_k,rv_v2_k,ri_k,grad_qv_k,...
    v_rows_k,iS_k,iS_ri_hat_k)
ri_k = ri_k+g_k;
g2_k = g_k.^2;
if ~isempty(v_k)
    v2_k = v_k.^2;
    % v^2*rv = t*qv*v^2 - v
    rv_v2_k = t*grad_qv_k.*v2_k-v_k;
    ri_k(v_rows_k,:) = ri_k(v_rows_k,:) + rv_v2_k;
end

iS_k = g2_k;
if ~isempty(v_k), iS_k(v_rows_k,:) = iS_k(v_rows_k,:) + v2_k; end
iS_k  = 1./iS_k;
iS_ri_hat_k = iS_k.*ri_k;
end

function [g2_ter,v2_ter,rv_v2_ter,ri_ter,iS_ter,iS_ri_hat_ter] = ...
    reduced_kkt_terminal_local(t,g_ter,v_ter,g2_ter,v2_ter,rv_v2_ter,...
    ri_ter,grad_qv_ter,iS_ter,iS_ri_hat_ter)
g2_ter = g_ter.^2;
ri_ter = ri_ter+g_ter;
if ~isempty(v_ter)
    v2_ter = v_ter.^2;
    % v^2*rv = t*qv*v^2 - v
    rv_v2_ter = t*grad_qv_ter.*v2_ter-v_ter;
    ri_ter = ri_ter + rv_v2_ter;
end

iS_ter = g2_ter;
if ~isempty(v_ter), iS_ter = iS_ter + v2_ter; end
iS_ter  = 1./iS_ter;
iS_ri_hat_ter = iS_ter.*ri_ter;
end
