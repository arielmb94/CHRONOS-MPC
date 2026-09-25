% UPDATE_MPC_STATE_CNSTR_LIMITS Update predicted-state bounds.
%
%   mpc = UPDATE_MPC_STATE_CNSTR_LIMITS(mpc, x_min, x_max) updates
%
%       x_min_k <= s_k <= x_max_k.
%
%   Use [] to leave either bound unchanged. The corresponding lower or
%   upper bound must first be enabled with INIT_MPC_STATE_CNSTR.
%
%   A bound may be a scalar, an nx-by-1 vector, or a time-varying nx-by-L
%   matrix, where L is the number of supplied horizon stages. If L < N, the
%   last supplied column is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     x_min   - Optional updated lower bound: scalar, nx-by-1, or nx-by-L.
%     x_max   - Optional updated upper bound: scalar, nx-by-1, or nx-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the upper bound:
%
%       mpc = update_mpc_state_cnstr_limits(mpc, [], x_max);
function mpc = update_mpc_state_cnstr_limits(mpc,min,max)

if ~isempty(min) && ~isempty(mpc.s_cnstr.min_limit)

    if isscalar(min)
        mpc.s_cnstr.min(:,:) = min;
    else
        mpc.s_cnstr.min(:,:) = fill_vec(mpc.s_cnstr.min, min, 1);
    end

end

if ~isempty(max) && ~isempty(mpc.s_cnstr.max_limit)
    
    if isscalar(max)
        mpc.s_cnstr.max(:,:) = max;
    else
        mpc.s_cnstr.max(:,:) = fill_vec(mpc.s_cnstr.max, max, 1);
    end

end

end
