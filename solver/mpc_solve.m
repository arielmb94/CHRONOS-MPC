% MPC_SOLVE Solve the current MPC problem.
%
%   [u0, mpc, iter] = MPC_SOLVE(mpc, s_prev, u_prev, r_in, xN_ref_in, ...
%       d_in, dz_in, dh_in) solves the current MPC problem.
%
%   Inputs:
%     mpc       - CHRONOS MPC structure with the persistent stage-local
%                 iterate.
%     s_prev    - nx column vector, last measured or estimated system state
%                 value.
%     u_prev    - nu column vector, control action applied to the system on
%                 the previous sampling time.
%     r_in      - Optional tracking reference. It can be an ny-by-1 column
%                 vector or an ny-by-L matrix, where L is the number of
%                 supplied horizon stages. A single column is reused for all
%                 stages. If L < N, the last supplied column is reused for
%                 the remaining stages; columns beyond the horizon are
%                 ignored. Pass [] when tracking is not used.
%     xN_ref_in - Optional nx column vector, reference for the terminal state
%                 xN of the prediction horizon. Required whenever the MPC
%                 problem contains terminal ingredients. If not used, pass
%                 an empty vector [].
%     d_in      - Optional known input for the dynamics and output model. It
%                 can be an nd-by-1 column vector or an nd-by-L matrix. A
%                 single column is reused for all stages. If L < N, the last
%                 supplied column is reused for the remaining stages; columns
%                 beyond the horizon are ignored. Pass [] when this input is
%                 not used.
%     dz_in     - Optional known input for the custom-cost signal z. It can
%                 be an ndz-by-1 column vector or an ndz-by-L matrix. A
%                 single column is reused for all stages. If L < N, the last
%                 supplied column is reused for the remaining stages; columns
%                 beyond the horizon are ignored. Pass [] when this input is
%                 not used.
%     dh_in     - Optional known input for the custom-constraint signal h. It
%                 can be an ndh-by-1 column vector or an ndh-by-L matrix. A
%                 single column is reused for all stages. If L < N, the last
%                 supplied column is reused for the remaining stages; columns
%                 beyond the horizon are ignored. Pass [] when this input is
%                 not used.
%
%   Output:
%     u0   - nu column vector, first step of the control action sequence
%            computed as the solution to the MPC problem.
%     mpc  - Updated CHRONOS MPC structure. Retain it and pass it to the
%            next call to MPC_SOLVE.
%     iter - Number of iterations required for the MPC optimization problem.
%
%   Example - solve with no optional runtime signals:
%
%       [u0, mpc, iter] = mpc_solve(mpc, s_prev, u_prev, [], [], [], [], []);
function [u0,mpc,iter] = mpc_solve(mpc,s_prev,u_prev,r_in,xN_ref_in,...
                                   d_in,dz_in,dh_in)

% handle input vector sizes
len_r_in = size(r_in,2);
if ~isempty(r_in)
    mpc.r(:,:) = fill_vec(mpc.r,r_in,1);
    if ~isempty(mpc.y_use_k0), mpc.r_0(:) = r_in(mpc.y_rows_k0,1); end
    if ~isempty(mpc.y_use_ter)
        ter_col = len_r_in;
        if ter_col > mpc.N, ter_col = mpc.N; end
        mpc.r_ter(:) = r_in(mpc.y_rows_ter,ter_col);
    end
end

if ~isempty(d_in)
    mpc.d(:,:) = fill_vec(mpc.d,d_in,1);
end

if ~isempty(dz_in)
    mpc.dz(:,:) = fill_vec(mpc.dz,dz_in,1);
end

if ~isempty(dh_in)
    mpc.dh(:,:) = fill_vec(mpc.dh,dh_in,1);
end

if ~isempty(mpc.ter_ingredients)
    if ~isempty(mpc.xN_ref_is_y) && isempty(xN_ref_in)
        mpc.xN_ref(:) = mpc.r_ter;
    else
        mpc.xN_ref(:) = xN_ref_in;
    end
end

% Recompute gradient/hessian if cost terms have been updated
if ~isempty(mpc.tracking_cost)
    if mpc.update_tracking
        mpc = update_tracking_cost(mpc,mpc.update_tracking);
    end
end
if ~isempty(mpc.quad_custom_cost)
    if mpc.update_customcost_quad
        mpc = update_custom_cost_quad(mpc,mpc.update_customcost_quad);
    end
end
if ~isempty(mpc.lin_custom_cost)
    if mpc.update_customcost_lin
        mpc = update_custom_cost_lin(mpc,mpc.update_customcost_lin);
    end
end
if mpc.recompute_cost_hess
    mpc = update_mpc_f0_hess(mpc,mpc.tracking_cost,...
        mpc.quad_control_cost,mpc.controlrate_cost,mpc.quad_custom_cost,...
        mpc.ter_ingredients);
end

% update dynamics equality RHS
mpc = update_mpc_beq(mpc,s_prev,mpc.A(:,:,1),mpc.dyn_use_d);

% Set Newton solver condition at start
continue_Newton = true;
iter = 0;

mpc = get_mpc_variables(mpc,mpc.has_du,mpc.tracking_cost,mpc.has_y_cnstr,...
                         mpc.has_h_cnstr,mpc.quad_custom_cost,mpc.lin_custom_cost,...
                         s_prev,u_prev);

lambda2 = 1;

while mpc.eps <= lambda2*0.5 && continue_Newton && iter < mpc.max_iter

    mpc = grad_f0_MPC(mpc,mpc.quad_control_cost,mpc.lin_control_cost,...
        mpc.controlrate_cost,mpc.tracking_cost,mpc.quad_custom_cost,...
        mpc.lin_custom_cost,mpc.ter_ingredients);

    mpc = equality_residuals(mpc,mpc.has_du);

    if ~isempty(mpc.g_0) || ~isempty(mpc.g_k) || ~isempty(mpc.g_ter)
        mpc = inequality_residuals(mpc,mpc.g_0,mpc.g_k,mpc.g_ter,...
            mpc.has_s_cnstr,mpc.has_u_cnstr,mpc.has_du_cnstr,...
            mpc.has_y_cnstr,mpc.has_h_cnstr);
    end

    mpc = reduced_KKT_elements(mpc,mpc.has_u_cnstr,mpc.has_du_cnstr,...
        mpc.has_s_cnstr,mpc.has_y_cnstr,mpc.has_h_cnstr,...
        mpc.g_0,mpc.g_k,mpc.g_ter);

    if ~isempty(mpc.has_du)
        mpc = riccati_KKT_du(mpc);
    else
        mpc = riccati_KKT(mpc);
    end

    if ~isempty(mpc.g_0) || ~isempty(mpc.g_k) || ~isempty(mpc.g_ter)
        mpc = recover_slacks(mpc,mpc.has_s_cnstr,mpc.has_u_cnstr,...
            mpc.has_du_cnstr,mpc.has_y_cnstr,mpc.has_h_cnstr);
    end

    lambda2 = get_lambda2(mpc);

    % Feasibility line search
    [mpc.g_0,mpc.g_k,mpc.g_ter,mpc.v_0,mpc.v_k,mpc.v_ter,mpc.u,mpc.s,...
        mpc.su,continue_Newton] = line_search_local(mpc.g_0,mpc.g_k,...
        mpc.g_ter,mpc.v_0,mpc.v_k,mpc.v_ter,mpc.u,mpc.s,mpc.su,...
        mpc.delta_g_0,mpc.delta_g_k,mpc.delta_g_ter,...
        mpc.delta_v_0,mpc.delta_v_k,mpc.delta_v_ter,mpc.delta_u,...
        mpc.delta_se,mpc.s_col,mpc.su_col,mpc.has_du,mpc.Beta,mpc.min_l,...
        continue_Newton);
    mpc = get_mpc_variables(mpc,mpc.has_du,mpc.tracking_cost,mpc.has_y_cnstr,...
                             mpc.has_h_cnstr,mpc.quad_custom_cost,mpc.lin_custom_cost,...
                             s_prev,u_prev);
    iter = iter+1;
end

u0 = mpc.u(:,1);

end

function [g_0,g_k,g_ter,v_0,v_k,v_ter,u,s,su,continue_Newton] = ...
    line_search_local(g_0,g_k,g_ter,v_0,v_k,v_ter,u,s,su,...
    delta_g_0,delta_g_k,delta_g_ter,delta_v_0,delta_v_k,delta_v_ter,...
    delta_u,delta_se,s_col,su_col,has_du,Beta,min_l,continue_Newton)
l = 1;
g_0_hat = g_0+l*delta_g_0;
g_k_hat = g_k+l*delta_g_k;
g_ter_hat = g_ter+l*delta_g_ter;
v_0_hat = v_0+l*delta_v_0;
v_k_hat = v_k+l*delta_v_k;
v_ter_hat = v_ter+l*delta_v_ter;

feas = all(g_0_hat(:)>0) && all(g_k_hat(:)>0) && all(g_ter_hat(:)>0) && ...
       all(v_0_hat(:)>0) && all(v_k_hat(:)>0) && all(v_ter_hat(:)>0);

if ~feas
    while ~feas
        l = l*Beta;

        g_0_hat = g_0+l*delta_g_0;
        g_k_hat = g_k+l*delta_g_k;
        g_ter_hat = g_ter+l*delta_g_ter;
        v_0_hat = v_0+l*delta_v_0;
        v_k_hat = v_k+l*delta_v_k;
        v_ter_hat = v_ter+l*delta_v_ter;

        feas = all(g_0_hat(:)>0) && all(g_k_hat(:)>0) && all(g_ter_hat(:)>0) && ...
               all(v_0_hat(:)>0) && all(v_k_hat(:)>0) && all(v_ter_hat(:)>0);
    end
    if l<min_l
        continue_Newton = false;
    end
end

g_0(:) = g_0_hat;
g_k(:,:) = g_k_hat;
g_ter(:) = g_ter_hat;
v_0(:) = v_0_hat;
v_k(:,:) = v_k_hat;
v_ter(:) = v_ter_hat;
u(:,:) = u+l*delta_u;
s(:,:) = s+l*delta_se(s_col,:);
if ~isempty(has_du)
    su(:,:) = su+l*delta_se(su_col,:);
end
end
