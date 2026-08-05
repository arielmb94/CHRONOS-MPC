% UPDATE_MPC_CONTROL_COST Update control-action cost weights.
%
%   mpc = UPDATE_MPC_CONTROL_COST(mpc, Ru, ru) replaces the quadratic
%   and/or linear weights in
%
%       J_control += 0.5*u_k'*Ru_k*u_k + ru_k'*u_k.
%
%   Use [] to leave either weight unchanged. The corresponding quadratic
%   or linear term must first be enabled with INIT_MPC_CONTROL_COST.
%
%   Ru may be nu-by-nu or nu-by-nu-by-L, and ru may be nu-by-1 or nu-by-L.
%   L is the number of supplied horizon stages. If L < N, the last supplied
%   stage is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Ru      - Optional updated quadratic weight, size nu-by-nu or
%               nu-by-nu-by-L.
%     ru      - Optional updated linear weight, size nu-by-1 or nu-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the linear weight:
%
%       mpc = update_mpc_Control_cost(mpc, [], ru);
function mpc = update_mpc_Control_cost(mpc,Ru,ru)

if ~isempty(Ru) && ~isempty(mpc.quad_control_cost)
    mpc.recompute_cost_hess = 1;
    mpc.Ru(:,:,:) = fill_mat(mpc.Ru, Ru, 1);
end

if ~isempty(ru) && ~isempty(mpc.lin_control_cost)
    mpc.ru(:,:) = fill_vec(mpc.ru, ru, 1);
end

end
