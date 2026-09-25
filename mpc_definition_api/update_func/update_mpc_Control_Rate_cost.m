% UPDATE_MPC_CONTROL_RATE_COST Update the quadratic control-rate weight.
%
%   mpc = UPDATE_MPC_CONTROL_RATE_COST(mpc, Rdu) replaces the weight in
%
%       J_rate += 0.5*Delta_u_k'*Rdu_k*Delta_u_k.
%
%   The control-rate cost must first be enabled with
%   INIT_MPC_CONTROLRATE_COST. Rdu may be nu-by-nu or nu-by-nu-by-L, where
%   L is the number of supplied horizon stages. If L < N, the last supplied
%   stage is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Rdu     - Updated control-rate weight, size nu-by-nu or
%               nu-by-nu-by-L. Each stage must be symmetric positive
%               semidefinite.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
function mpc = update_mpc_Control_Rate_cost(mpc,Rdu)

if ~isempty(mpc.controlrate_cost)
    mpc.recompute_cost_hess = 1;
    mpc.Rdu(:,:,:) = fill_mat(mpc.Rdu, Rdu, 1);
end

end
