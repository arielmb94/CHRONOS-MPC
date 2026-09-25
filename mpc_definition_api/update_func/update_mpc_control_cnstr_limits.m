% UPDATE_MPC_CONTROL_CNSTR_LIMITS Update control-action bounds.
%
%   mpc = UPDATE_MPC_CONTROL_CNSTR_LIMITS(mpc, u_min, u_max) updates
%
%       u_min_k <= u_k <= u_max_k.
%
%   Use [] to leave either bound unchanged. The corresponding lower or
%   upper bound must first be enabled with INIT_MPC_CONTROL_CNSTR.
%
%   A bound may be a scalar, an nu-by-1 vector, or a time-varying nu-by-L
%   matrix, where L is the number of supplied horizon stages. If L < N, the
%   last supplied column is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     u_min   - Optional updated lower bound: scalar, nu-by-1, or nu-by-L.
%     u_max   - Optional updated upper bound: scalar, nu-by-1, or nu-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update both bounds:
%
%       mpc = update_mpc_control_cnstr_limits(mpc, u_min, u_max);
function mpc = update_mpc_control_cnstr_limits(mpc,min,max)

if ~isempty(min) && ~isempty(mpc.u_cnstr.min_limit)

    if isscalar(min)
        mpc.u_cnstr.min(:,:) = min;
    else
        mpc.u_cnstr.min(:,:) = fill_vec(mpc.u_cnstr.min, min, 1);
    end
    
end

if ~isempty(max) && ~isempty(mpc.u_cnstr.max_limit)

    if isscalar(max)
        mpc.u_cnstr.max(:,:) = max;
    else
        mpc.u_cnstr.max(:,:) = fill_vec(mpc.u_cnstr.max, max, 1);
    end

end

end
