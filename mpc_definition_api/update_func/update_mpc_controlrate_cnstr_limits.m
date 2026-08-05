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
function mpc = update_mpc_controlrate_cnstr_limits(mpc,min,max)

if ~isempty(min)

    if isscalar(min)
        mpc.du_cnstr.min(:,:) = min;
    elseif size(min,2) < mpc.N
        mpc.du_cnstr.min(:,:) = fill_vec(mpc.du_cnstr.min, min, 1);
    else
        mpc.du_cnstr.min(:,:) = min(:,1:mpc.N);
    end
    
    mpc.bi_0(mpc.du_cnstr.min_ineqRow_0) = -mpc.du_cnstr.min(:,1);
    mpc.bi_k(mpc.du_cnstr.min_ineqRow_k,:) = -mpc.du_cnstr.min(:,2:mpc.N);
end

if ~isempty(max)

    if isscalar(max)
        mpc.du_cnstr.max(:,:) = max;
    elseif size(max,2) < mpc.N
        mpc.du_cnstr.max(:,:) = fill_vec(mpc.du_cnstr.max, max, 1);
    else
        mpc.du_cnstr.max(:,:) = max(:,1:mpc.N);
    end

    mpc.bi_0(mpc.du_cnstr.max_ineqRow_0) = mpc.du_cnstr.max(:,1);
    mpc.bi_k(mpc.du_cnstr.max_ineqRow_k,:) = mpc.du_cnstr.max(:,2:mpc.N);
end

end