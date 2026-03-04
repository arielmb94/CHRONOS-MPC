% INIT_INITIAL_GUESS Generates a feasible warm-start primal vector for the MPC solver.
%
%   x0 = INIT_INITIAL_GUESS(mpc, s_prev, u_prev) calculates a strictly feasible 
%   initial guess for the primal optimization vector using a system rollout.
%   It mathematically guarantees that the initial guess respects input limits, 
%   rate limits, and hard state constraints, preventing solver crashes.
%
%   x0 = INIT_INITIAL_GUESS(mpc, s_prev, u_prev, x_ref, d_in, dh_in) allows the 
%   inclusion of reference trajectories and measured disturbances.
%
%   HOW IT WORKS:
%   The function simulates the system dynamics across the prediction horizon:
%     1. Control Strategy: If an LQR gain (mpc.K) and reference (x_ref) exist, 
%        it uses feedback (u = K*(x_ref - x)). Otherwise, it holds u_prev constant.
%     2. Input Clipping: Control actions are strictly bounded by absolute (mpc.u_cnstr)
%        and rate (mpc.du_cnstr) limits.
%     3. State Clamping: If hard state constraints (mpc.s_cnstr) are defined (without 
%        slacks), the state is clamped to remain strictly feasible, intentionally 
%        breaking the equality constraint to prevent an Interior Point method crash.
%     4. Soft Constraints: States with enabled slacks are allowed to violate bounds. 
%        The slack variables are then automatically sized to absorb the violation.
%
%   INPUTS:
%       mpc    - CHRONOS MPC structure.
%       s_prev - [nx x 1] Current measured state vector.
%       u_prev - [nu x 1] Last applied control input.
%       x_ref  - [nx x 1] (Optional) State reference. Defaults to [].
%       d_in   - [nd x 1] (Optional) Measured disturbances.
%       dh_in  - [ndh x 1] (Optional) Measured disturbance vector for custom constraints.
%
%   OUTPUTS:
%       x0     - Primal vector [u_0; x_1; u_1; ... x_Nc; ... x_N; v] 
function x0 = init_initial_guess(mpc,s_prev,u_prev,x_ref,d_in,dh_in)
arguments
    mpc
    s_prev
    u_prev
    x_ref = []
    d_in = []
    dh_in = []
end
    
    x0 = zeros(mpc.Nx+mpc.Nu+mpc.Nv,1);

    % handle input vectors size
    if ~isempty(d_in) && length(d_in)< mpc.Nd
        d = fill_vec(d_in,mpc.nd,mpc.Nd,1);
    else
        d = d_in;
    end

    if ~isempty(dh_in) && length(dh_in)< mpc.Nd
        dh = fill_vec(dh_in,mpc.ndh,mpc.Ndh,1);
    else
        dh = dh_in;
    end
    
    x_k = s_prev;
    u_k_prev = u_prev;
    
    for k = 1 : mpc.N_ctr_hor
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
        if ~isempty(mpc.du_cnstr)
            % Check minimum rate limit
            if ~isempty(mpc.du_cnstr.min)
                du_min_strict = mpc.du_cnstr.min + mpc.slack_epsilon;
                u_k = max(u_k_prev + du_min_strict, u_k);
            end
            % Check maximum rate limit
            if ~isempty(mpc.du_cnstr.max)
                du_max_strict = mpc.du_cnstr.max - mpc.slack_epsilon;
                u_k = min(u_k_prev + du_max_strict, u_k);
            end
        end
        
        % 3. Clip for Absolute Constraints (u)
        if ~isempty(mpc.u_cnstr)
            % Check minimum absolute limit
            if ~isempty(mpc.u_cnstr.min)
                u_min_strict = mpc.u_cnstr.min + mpc.slack_epsilon;
                u_k = max(u_min_strict, u_k);
            end
            % Check maximum absolute limit
            if ~isempty(mpc.u_cnstr.max)
                u_max_strict = mpc.u_cnstr.max - mpc.slack_epsilon;
                u_k = min(u_max_strict, u_k);
            end
        end
        
        % 4. Propagate Dynamics
        % x_{k+1} = A*x_k + B*u_k + D*d_k
        % Assuming d_seq is [nd x N]
        x_next = mpc.A * x_k + mpc.B * u_k;
        if ~isempty(mpc.Bd) && ~isempty(d_in)
            x_next = x_next + mpc.Bd * d_in;
        end
        % clamp x
        if ~isempty(mpc.s_cnstr)
            % Check minimum state limits
            if ~isempty(mpc.s_cnstr.min)
                % Only clamp if this constraint DOES NOT have a slack variable
                if ~mpc.s_cnstr.min_slack_nv 
                    s_min_strict = mpc.s_cnstr.min + mpc.slack_epsilon;
                    x_next = max(s_min_strict, x_next);
                end
            end
            
            % Check maximum state limits
            if ~isempty(mpc.s_cnstr.max)
                % Only clamp if this constraint DOES NOT have a slack variable
                if ~mpc.s_cnstr.max_slack_nv
                    s_max_strict = mpc.s_cnstr.max - mpc.slack_epsilon;
                    x_next = min(s_max_strict, x_next);
                end
            end
        end
        
        % 5. Map to Primal Optimization Vector
        % Indexing math for [u0; x1; u1; x2; ...]
        idx_start = (k - 1) * (mpc.nu + mpc.nx);
        
        % Insert u_k
        x0(idx_start + 1 : idx_start + mpc.nu) = u_k;
        
        % Insert x_{next}
        x0(idx_start + mpc.nu + 1 : idx_start + mpc.nu + mpc.nx) = x_next;
        
        % 6. Prepare for next step
        x_k = x_next;
        u_k_prev = u_k;
    end
    
    if mpc.N_ctr_hor < mpc.N
        u_N_ctr_hor = u_k;
        idx_start_N_ctr_hor = (mpc.N_ctr_hor-1)*(mpc.nu + mpc.nx) + mpc.nu;
        for k = mpc.N_ctr_hor + 1 : mpc.N
            % 4. Propagate Dynamics
            % x_{k+1} = A*x_k + B*u_k + D*d_k
            x_next = mpc.A * x_k + mpc.B * u_N_ctr_hor;
            if ~isempty(mpc.Bd) && ~isempty(d_in)
                x_next = x_next + mpc.Bd * d_in;
            end
            % clamp x
            if ~isempty(mpc.s_cnstr)
                % Check minimum state limits
                if ~isempty(mpc.s_cnstr.min)
                    % Only clamp if this constraint DOES NOT have a slack variable
                    if ~mpc.s_cnstr.min_slack_nv
                        s_min_strict = mpc.s_cnstr.min + mpc.slack_epsilon;
                        x_next = max(s_min_strict, x_next);
                    end
                end

                % Check maximum state limits
                if ~isempty(mpc.s_cnstr.max)
                    % Only clamp if this constraint DOES NOT have a slack variable
                    if ~mpc.s_cnstr.max_slack_nv
                        s_max_strict = mpc.s_cnstr.max - mpc.slack_epsilon;
                        x_next = min(s_max_strict, x_next);
                    end
                end
            end

            % 5. Map to Primal Optimization Vector
            % Indexing math for [u0; x1; u1; x2; ...]
            idx_start = idx_start_N_ctr_hor + (k - mpc.N_ctr_hor) * mpc.nx;

            % Insert x_{next}
            x0(idx_start + 1 : idx_start + mpc.nx) = x_next;

            % 6. Prepare for next step
            x_k = x_next;
        end
    end

    % set x_ref to s_prev if user does not pass x_ref but there are
    % terminal constraints
    if isempty(x_ref) && mpc.ter_constraint
        x_ref = s_prev;
    end
    
    if mpc.Nv
        x0(mpc.Nx+mpc.Nu+1:mpc.Nx+mpc.Nu+mpc.Nv) = mpc.slack_epsilon*2;
        mpc = get_mpc_variables(mpc,x0,s_prev,u_prev,[],d,dh,[]);
        [~,x0,feas] = mpc_slack_update(mpc,x0,x_ref);
        if ~feas
            msg = ['Cannot find feasible initial guess due to hard constraints, ' ...
                'enable slack variables on output/custom constraints'];
            warning(msg);
        end
    else
        mpc = get_mpc_variables(mpc,x0,s_prev,u_prev,[],d,dh,[]);
        [~,feas] = check_mpc_feasibility(mpc,x_ref);
        if ~feas
            msg = ['Cannot find feasible initial guess due to hard constraints, ' ...
                'enable slack variables on output/custom constraints'];
            warning(msg);
        end
    end
end