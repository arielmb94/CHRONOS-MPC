% UPDATE_MPC_SLACK_COST Update soft-constraint violation penalties.
%
%   mpc = UPDATE_MPC_SLACK_COST(mpc, cnstr, qv_min, qv_max) replaces the
%   lower- and/or upper-bound violation penalties for an initialized soft
%   constraint. Pass one of these structures as cnstr:
%
%       mpc.s_cnstr    predicted-state constraints
%       mpc.y_cnstr    output constraints
%       mpc.h_cnstr    custom-signal constraints
%
%   Use [] to leave either penalty unchanged. The corresponding lower or
%   upper bound must already exist. Larger qv values make violations more
%   costly, while the constraint remains soft and may still be violated
%   when it cannot be satisfied.
%
%   A penalty may be a scalar, an n-by-1 vector, or a time-varying n-by-L
%   matrix, where n is nx, ny, or nh for the selected constraint and L is
%   the number of supplied horizon stages. If L < N, the last supplied
%   column is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     cnstr   - Initialized mpc.s_cnstr, mpc.y_cnstr, or mpc.h_cnstr.
%     qv_min  - Optional lower-bound penalty: scalar, n-by-1, or n-by-L.
%     qv_max  - Optional upper-bound penalty: scalar, n-by-1, or n-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only state lower-bound penalties:
%
%       mpc = update_mpc_slack_cost(mpc, mpc.s_cnstr, qv_min, []);
function mpc = update_mpc_slack_cost(mpc, cnstr, qv_min, qv_max)
% Note: The struct fields cnstr.qv_min and cnstr.qv_max are intentionally 
% NOT updated here to maintain a clean single-output API (returning only mpc). 
% The underlying solver only requires the updated global mpc.grad_qv vectors.

if ~isempty(qv_min) && ~isempty(cnstr) && ~isempty(cnstr.min_limit)
    if isscalar(qv_min)

        if ~isempty(cnstr.use_k0)
            mpc.grad_qv_0(cnstr.min_row_v_0) = qv_min; 
        end
        
        mpc.grad_qv_k(cnstr.min_row_v_k, :) = qv_min;
        
        if ~isempty(cnstr.use_ter)
            mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min; 
        end
    else
        % Initial stage (k = 0)
        if ~isempty(cnstr.use_k0)
            mpc.grad_qv_0(cnstr.min_row_v_0) = qv_min(cnstr.rows_k0, 1); 
        end
        
        % Intermediate stages (k = 1 ... N-1)
        mpc.grad_qv_k(cnstr.min_row_v_k, :) = fill_vec(mpc.grad_qv_k(cnstr.min_row_v_k, :), qv_min, 1);
        
        % Terminal stage (k = N)
        if ~isempty(cnstr.use_ter)
            ter_col = size(qv_min, 2);
            if ter_col > mpc.N, ter_col = mpc.N; end
            mpc.grad_qv_ter(cnstr.min_ineqRow_ter) = qv_min(cnstr.rows_ter, ter_col);
        end
    end
end

if ~isempty(qv_max) && ~isempty(cnstr) && ~isempty(cnstr.max_limit)
    if isscalar(qv_max)
        % Fully vectorized scalar assignment
        if ~isempty(cnstr.use_k0)
            mpc.grad_qv_0(cnstr.max_row_v_0) = qv_max; 
        end
        
        mpc.grad_qv_k(cnstr.max_row_v_k, :) = qv_max;
        
        if ~isempty(cnstr.use_ter)
            mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max; 
        end
    else
        % Initial stage (k = 0)
        if ~isempty(cnstr.use_k0)
            mpc.grad_qv_0(cnstr.max_row_v_0) = qv_max(cnstr.rows_k0, 1); 
        end
        
        % Intermediate stages (k = 1 ... N-1)
        mpc.grad_qv_k(cnstr.max_row_v_k, :) = fill_vec(mpc.grad_qv_k(cnstr.max_row_v_k, :), qv_max, 1);
        
        % Terminal stage (k = N)
        if ~isempty(cnstr.use_ter)
            ter_col = size(qv_max, 2);
            if ter_col > mpc.N, ter_col = mpc.N; end
            mpc.grad_qv_ter(cnstr.max_ineqRow_ter) = qv_max(cnstr.rows_ter, ter_col);
        end
    end
end
end
