% BUILD_CHRONOS_MPC This function must always be called at the end of
% the initialization process. Initializes a feasible stage-local iterate for
% the MPC solver and initializes costs and constraint structures based on
% initialization data. 
%
%   mpc = BUILD_CHRONOS_MPC(mpc, s_prev, u_prev) calculates a strictly feasible
%   initial iterate using a system rollout.
%
%   mpc = BUILD_CHRONOS_MPC(mpc, s_prev, u_prev, d_in, dh_in, x_ref) allows the
%   inclusion of reference trajectories and measured disturbances.
%
%   Inputs:
%     mpc    - CHRONOS MPC structure.
%     s_prev - [nx x 1] Current measured state vector.
%     u_prev - [nu x 1] Last applied control input.
%     d_in   - Optional known input for the dynamics and output model. It
%              can be an nd-by-1 column vector or an nd-by-L matrix, where
%              L is the number of supplied horizon stages. A single column
%              is reused across the horizon. If L < N, the last supplied
%              column is reused for the remaining stages; columns beyond
%              the horizon are ignored. Pass [] when this input is not used.
%     dh_in  - Optional known input for the custom-constraint signal. It can
%              be an ndh-by-1 column vector or an ndh-by-L matrix, where L
%              is the number of supplied horizon stages. A single column is
%              reused across the horizon. If L < N, the last supplied
%              column is reused for the remaining stages; columns beyond
%              the horizon are ignored. Pass [] when this input is not used.
%     x_ref  - [nx x 1] (Optional) State reference. Defaults to [].
%
%   Output:
%     mpc    - MPC structure with initialized stage-local primal fields.
%
%   Example - build with no optional runtime signals:
%
%       mpc = build_chronos_mpc(mpc, s_prev, u_prev, [], [], []);
function mpc = build_chronos_mpc(mpc,s_prev,u_prev,d_in,dh_in,x_ref)
arguments
    mpc
    s_prev
    u_prev
    d_in = []
    dh_in = [] 
    x_ref = []
end

    % allocate stage-local optimization variables and solver workspace
    mpc = build_optimization_variables(mpc);
    
    % init equality constraints
    mpc = genEqualities(mpc);

    % assign constraints rows
    mpc = assign_inequalities(mpc);

    % init soft slacks cost
    mpc = set_soft_cost_qv(mpc);

    % init costs
    mpc = init_costs(mpc);

    % preallocate fixed-size Riccati workspace for the online Newton solve
    mpc = preallocate_riccati(mpc);

    % initialize stage-local primal variables
    if ~isempty(d_in)
        mpc.d(:,:) = fill_vec(mpc.d,d_in,1);
    end

    mpc = rollstates(mpc,s_prev,u_prev,x_ref,mpc.d);

    if ~isempty(dh_in)
        mpc.dh(:,:) = fill_vec(mpc.dh,dh_in,1);
    end

    mpc = get_mpc_variables(mpc,mpc.has_du,mpc.tracking_cost,mpc.has_y_cnstr,...
                            mpc.has_h_cnstr,mpc.quad_custom_cost,mpc.lin_custom_cost,...
                            s_prev,u_prev);
    mpc = initialize_inequality_slacks(mpc);
    
end

function mpc = initialize_inequality_slacks(mpc)

if isempty(mpc.g_0) && isempty(mpc.g_k) && isempty(mpc.g_ter)
    return
end

% Form raw bound residuals through the production family and row-map path.
mpc.g_0(:) = 0;
mpc.g_k(:,:) = 0;
mpc.g_ter(:) = 0;
mpc.v_0(:) = 0;
mpc.v_k(:,:) = 0;
mpc.v_ter(:) = 0;
mpc = inequality_residuals(mpc,mpc.g_0,mpc.g_k,mpc.g_ter,...
    mpc.has_s_cnstr,mpc.has_u_cnstr,mpc.has_du_cnstr,...
    mpc.has_y_cnstr,mpc.has_h_cnstr);

[mpc.g_0,mpc.v_0] = initialize_inequality_slacks_local(...
    mpc.g_0,mpc.v_0,mpc.ri_0,mpc.v_rows_0,mpc.slack_epsilon);
[mpc.g_k,mpc.v_k] = initialize_inequality_slacks_local(...
    mpc.g_k,mpc.v_k,mpc.ri_k,mpc.v_rows_k,mpc.slack_epsilon);

mpc.g_ter(:) = max(-mpc.ri_ter,mpc.slack_epsilon);
if ~isempty(mpc.v_ter)
    mpc.g_ter(:) = mpc.slack_epsilon + max(-mpc.ri_ter,0);
    mpc.v_ter(:) = mpc.slack_epsilon + max(mpc.ri_ter,0);
end

mpc = inequality_residuals(mpc,mpc.g_0,mpc.g_k,mpc.g_ter,...
    mpc.has_s_cnstr,mpc.has_u_cnstr,mpc.has_du_cnstr,...
    mpc.has_y_cnstr,mpc.has_h_cnstr);

end

function [g,v] = initialize_inequality_slacks_local(g,v,ri,v_rows,epsilon)

g(:,:) = max(-ri,epsilon);
if ~isempty(v)
    ri_soft = ri(v_rows,:);
    g(v_rows,:) = epsilon + max(-ri_soft,0);
    v(:,:) = epsilon + max(ri_soft,0);
end

end

function mpc = rollstates(mpc,s_prev,u_prev,x_ref,d)

x_k = s_prev;
u_k_prev = u_prev;
 
for k = 1 : mpc.N
    % Compute Raw Control Input
    if ~isempty(mpc.K) && ~isempty(x_ref)
        % Option 2: Terminal ingredients exist and user pass x_ref
        u_raw = mpc.K * (x_ref - x_k);
    else
        % Option 1: Constant previous input
        u_raw = u_k_prev;
    end

    u_k = u_raw;

% 2. Clip for Rate Constraints (Delta u)
if ~isempty(mpc.has_du_cnstr)
    % Check minimum rate limit
    if ~isempty(mpc.du_cnstr.min_limit)
        du_min_strict = mpc.du_cnstr.min(:,k) + mpc.slack_epsilon;
        u_k = max(u_k_prev + du_min_strict, u_k);
    end
    % Check maximum rate limit
    if ~isempty(mpc.du_cnstr.max_limit)
        du_max_strict = mpc.du_cnstr.max(:,k) - mpc.slack_epsilon;
        u_k = min(u_k_prev + du_max_strict, u_k);
    end
end

% 3. Clip for Absolute Constraints (u)
if ~isempty(mpc.has_u_cnstr)
    % Check minimum absolute limit
    if ~isempty(mpc.u_cnstr.min_limit)
        u_min_strict = mpc.u_cnstr.min(:,k) + mpc.slack_epsilon;
        u_k = max(u_min_strict, u_k);
    end
    % Check maximum absolute limit
    if ~isempty(mpc.u_cnstr.max_limit)
        u_max_strict = mpc.u_cnstr.max(:,k) - mpc.slack_epsilon;
        u_k = min(u_max_strict, u_k);
    end
end

    % 4. Propagate Dynamics
    % x_{k+1} = A*x_k + B*u_k + D*d_k
    x_next = mpc.A(:,:,k) * x_k + mpc.B(:,:,k) * u_k;
    if ~isempty(mpc.dyn_use_d) && any(d(:))
        x_next = x_next + mpc.Bd(:,:,k) * d(:,k);
    end
% clamp x: for safety net in case we are dealing with unstable
% system
if ~isempty(mpc.has_s_cnstr)
    % Check minimum state limits
    if ~isempty(mpc.s_cnstr.min_limit)
        s_min_strict = mpc.s_cnstr.min(:,k) + mpc.slack_epsilon;
            x_next = max(s_min_strict, x_next);
    end

    % Check maximum state limits
    if ~isempty(mpc.s_cnstr.max_limit)
        s_max_strict = mpc.s_cnstr.max(:,k) - mpc.slack_epsilon;
            x_next = min(s_max_strict, x_next);

    end
end

    % 5. Store the stage-local iterate
    mpc.u(:,k) = u_k;
    if ~isempty(mpc.has_du)
        mpc.su(:,k) = u_k_prev;
    end
    mpc.s(:,k) = x_next;

    % 6. Prepare for next step
    x_k = x_next;
    u_k_prev = u_k;
end

end

function mpc = preallocate_riccati(mpc)
    
    n_interior = mpc.N-1;
    mpc.Q_hat = zeros(mpc.nse,mpc.nse,mpc.N);
    mpc.R_hat = zeros(mpc.nu,mpc.nu,n_interior);
    mpc.Y_hat = zeros(mpc.nu,mpc.nse,n_interior);
    mpc.K_ric = zeros(mpc.nu,mpc.nse,n_interior);
    mpc.d_ric = zeros(mpc.nu,n_interior);
    mpc.rs_hat = zeros(mpc.nse,mpc.N);
    mpc.rp_hat = zeros(mpc.nse,n_interior);
    mpc.ru_hat = zeros(mpc.nu,n_interior);

    mpc.R_hat_0 = zeros(mpc.nu,mpc.nu);
    mpc.d_ric_0 = zeros(mpc.nu,1);
    mpc.rp_hat_0 = zeros(mpc.nse,1);
    mpc.ru_hat_0 = zeros(mpc.nu,1);

    mpc.QA_ric = zeros(mpc.nx,mpc.nx);
    if ~isempty(mpc.has_du)
        mpc.QB_ric = [];
        mpc.G_ric = zeros(mpc.nu,mpc.nx);
    else
        mpc.QB_ric = zeros(mpc.nx,mpc.nu);
        mpc.G_ric = [];
    end
    mpc.solve_rhs_ric = zeros(mpc.nu,mpc.nse+1);
    mpc.solve_result_ric = zeros(mpc.nu,mpc.nse+1);
end
