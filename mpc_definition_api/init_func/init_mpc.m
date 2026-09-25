% INIT_MPC Create the CHRONOS MPC structure.
%
%   mpc = INIT_MPC(N) creates an empty MPC structure with prediction
%   horizon N.
%
%   mpc = INIT_MPC() uses the default horizon N = 10.
%
%   Call this function first. Then define the dynamics, optionally replace
%   the default tracking output, add costs and constraints, and finish with
%   BUILD_CHRONOS_MPC.
%
%   Input:
%     N       - Prediction horizon. Default: 10.
%
%   Output:
%     mpc     - New CHRONOS MPC structure, ready for model, cost, and
%               constraint initialization.
%
%   Example:
%
%       mpc = init_mpc(20);
%       mpc = init_mpc_dynamics(mpc, A, B, []);
%       mpc = init_mpc_output(mpc, C, [], []);
function mpc = init_mpc(N)
arguments
    N = 10
end

mpc.N = N;

mpc.Qe = [];
mpc.Qe_0 = [];
mpc.Qe_ter = [];
mpc.Rdu = [];
mpc.Ru = [];
mpc.ru = [];
mpc.A = [];
mpc.B = [];
mpc.Bd = [];
mpc.C = [];
mpc.D = [];
mpc.Dd = [];
mpc.C_0 = [];
mpc.D_0 = [];
mpc.Dd_0 = [];
mpc.C_ter = [];
mpc.y_rows_k0 = [];
mpc.y_rows_ter = [];
mpc.Qz = [];
mpc.Qz_0 = [];
mpc.Qz_ter = [];
mpc.qz = [];
mpc.qz_0 = [];
mpc.qz_ter = [];
mpc.Cz = [];
mpc.Dz = [];
mpc.Ddz = [];
mpc.Dsuz = [];
mpc.Cz_0 = [];
mpc.Dz_0 = [];
mpc.Dsuz_0 = [];
mpc.Ddz_0 = [];
mpc.Cz_ter = [];
mpc.z_rows_k0 = [];
mpc.z_rows_ter = [];
mpc.Ch = [];
mpc.Dh = [];
mpc.Dsuh = [];
mpc.Ddh = [];
mpc.Dh_0 = [];
mpc.Dsuh_0 = [];
mpc.Ddh_0 = [];
mpc.Ch_ter = [];
mpc.r_0 = [];
mpc.r_ter = [];
mpc.nx = 0;
mpc.nu = 0;
mpc.nd = 0;
mpc.ny = 0;
mpc.ny_0 = 0;
mpc.ny_ter = 0;
mpc.ndz = 0;
mpc.nz_0 = 0;
mpc.nz = 0;
mpc.nz_ter = 0;
mpc.ndh = 0;
mpc.nh_0 = 0;
mpc.nh = 0;
mpc.nh_ter = 0;
mpc.dyn_use_d = [];
mpc.y_use_s = [];
mpc.y_use_u = [];
mpc.y_use_d = [];
mpc.y_use_k0 = [];
mpc.y_use_ter = [];
mpc.z_use_s = [];
mpc.z_use_u = [];
mpc.z_use_su = [];
mpc.z_use_d = [];
mpc.z_use_k0 = [];
mpc.z_use_ter = [];

mpc.s = [];
mpc.s_ter = [];
mpc.u = [];
mpc.su = [];
mpc.du = [];
mpc.r = [];
mpc.y = [];
mpc.err = [];
mpc.err_0 = [];
mpc.err_ter = [];
mpc.d = [];
mpc.z = [];
mpc.z_0 = [];
mpc.z_ter = [];
mpc.dz = [];
mpc.h = [];
mpc.dh = [];
mpc.g_0 = [];
mpc.g_k = [];
mpc.g_ter = [];
mpc.v_0 = [];
mpc.v_k = [];
mpc.v_ter = [];
mpc.v_rows_0 = [];
mpc.v_rows_k = [];
mpc.xN_ref = [];
mpc.y_0 = [];
mpc.y_ter = [];
mpc.h_0 = [];
mpc.h_ter = [];
mpc.beq_0 = [];
mpc.beq_k = [];
mpc.rp_0 = [];
mpc.rp_k = [];
mpc.ri_0 = [];
mpc.ri_k = [];
mpc.ri_ter = [];
mpc.delta_u = [];
mpc.delta_se = [];

mpc.tracking_cost = [];
mpc.quad_control_cost = [];
mpc.lin_control_cost = [];
mpc.controlrate_cost = [];
mpc.quad_custom_cost = [];
mpc.lin_custom_cost = [];
mpc.s_col = [];
mpc.su_col = [];
mpc.u_col = [];
mpc.se_col = [];

mpc.gradErr_Qe_0 = [];
mpc.gradErr_Qe_k = [];
mpc.gradErr_Qe_ter = [];
mpc.gradz_Qz_0 = [];
mpc.gradz_Qz_k = [];
mpc.gradz_Qz_ter = [];
mpc.gradz_qz_0 = [];
mpc.gradz_qz_k = [];
mpc.gradz_qz_ter = [];
mpc.gradRateCtrl_Rdu_k = [];

mpc.H_f0_0 = [];
mpc.H_f0_k = [];
mpc.H_f0_ter = [];
mpc.R_E_0 = [];
mpc.Q_E = [];
mpc.R_E = [];
mpc.Y_E = [];
mpc.Q_E_ter = [];
mpc.R_Z_0 = [];
mpc.Q_Z = [];
mpc.R_Z = [];
mpc.Y_Z = [];
mpc.Q_Z_ter = [];
mpc.Q_k = [];
mpc.Q_ter = [];
mpc.R_0 = [];
mpc.R_k = [];
mpc.Y_k = [];
mpc.grad_qv_0 = [];
mpc.grad_qv_k = [];
mpc.grad_qv_ter = [];
mpc.g2_0 = [];
mpc.g2_k = [];
mpc.g2_ter = [];
mpc.v2_0 = [];
mpc.v2_k = [];
mpc.v2_ter = [];
mpc.rv_v2_0 = [];
mpc.rv_v2_k = [];
mpc.rv_v2_ter = [];
mpc.iS_0 = [];
mpc.iS_k = [];
mpc.iS_ter = [];
mpc.iS_ri_hat_0 = [];
mpc.iS_ri_hat_k = [];
mpc.iS_ri_hat_ter = [];
mpc.ru_hat_0 = [];
mpc.delta_g_0 = [];
mpc.delta_g_k = [];
mpc.delta_g_ter = [];
mpc.delta_v_0 = [];
mpc.delta_v_k = [];
mpc.delta_v_ter = [];

mpc.grad_f0_0 = [];
mpc.grad_f0_k = [];
mpc.grad_f0_ter = [];
mpc.ru_0 = [];
mpc.rse_k = [];
mpc.ru_k = [];
mpc.rse_ter = [];

mpc.update_tracking = false;
mpc.update_customcost_quad = false;
mpc.update_customcost_lin = false;
mpc.recompute_cost_hess = false;

mpc.s = [];
mpc.s_ter = [];
mpc.u = [];
mpc.du = [];
mpc.y = [];
mpc.h = [];
mpc.z = [];
mpc.n = 0;
mpc.nvar = 0;
mpc.t = 50;
mpc.Beta = 0.75;
mpc.min_l = 1e-6;
mpc.eps = 1e-4;
mpc.max_iter = 10;
mpc.ter_ingredients = [];
mpc.xN_ref_is_y = [];
mpc.P = [];
mpc.P2 = [];
mpc.K = [];
mpc.warm_starting = false;

mpc.has_s_cnstr = [];
mpc.has_u_cnstr = [];
mpc.has_du = [];
mpc.has_du_cnstr = [];
mpc.has_y_cnstr = [];
mpc.has_h_cnstr = [];

mpc.s_cnstr = [];
mpc.u_cnstr = [];
mpc.du_cnstr = [];
mpc.y_cnstr = [];
mpc.h_cnstr = [];

mpc.ng_k = [0 0 0]; % inequalites per horizon step
mpc.nv_k = [0 0 0]; % soft inequalites per horizon step

mpc.qv = 50; % Slack variable penalty
mpc.slack_epsilon = 1e-3;
mpc.eps_thknv = 1e-6;
mpc.H_CustomCost_0 = [];

end
