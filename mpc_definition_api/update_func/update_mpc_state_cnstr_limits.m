<<<<<<< HEAD
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
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_u_cnstr(mpc,u_min,u_max)
%
% Modifies the constraints limits on the control action
%
% In:
%   - mpc: CHRONOS mpc structure
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_state_cnstr_limits(mpc,min,max)

if ~isempty(min)

    if isscalar(min)
        mpc.s_cnstr.min(:,:) = min;
    elseif size(min,2) < mpc.N
        mpc.s_cnstr.min(:,:) = fill_vec(mpc.s_cnstr.min, min, 1);
    else
        mpc.s_cnstr.min(:,:) = min(:,1:mpc.N);
    end

    mpc.bi_k(mpc.s_cnstr.min_ineqRow_k,:) = -mpc.s_cnstr.min(:,1:mpc.N-1);
    mpc.bi_ter(mpc.s_cnstr.min_ineqRow_ter) = -mpc.s_cnstr.min(:,mpc.N);
end

if ~isempty(max)
    
    if isscalar(max)
        mpc.s_cnstr.max(:,:) = max;
    elseif size(max,2) < mpc.N
        mpc.s_cnstr.max(:,:) = fill_vec(mpc.s_cnstr.max, max, 1);
    else
        mpc.s_cnstr.max(:,:) = max(:,1:mpc.N);
    end

    mpc.bi_k(mpc.s_cnstr.max_ineqRow_k,:) = mpc.s_cnstr.max(:,1:mpc.N-1);
    mpc.bi_ter(mpc.s_cnstr.max_ineqRow_ter) = mpc.s_cnstr.max(:,mpc.N);
end

end
>>>>>>> 9d5c3b4 (time varying constraint limits and slack costs)
