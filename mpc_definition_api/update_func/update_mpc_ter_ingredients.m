% UPDATE_MPC_TER_INGREDIENTS Update the terminal-state cost weight.
%
%   mpc = UPDATE_MPC_TER_INGREDIENTS(mpc, P) replaces the weight in
%
%       J_terminal = (xN_ref-s_N)'*P*(xN_ref-s_N).
%
%   This function accepts a new P directly; it does not rerun the 
%   terminal ingredient calculation.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     P       - Updated terminal weight, size nx-by-nx, symmetric positive
%               definite.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
function mpc = update_mpc_ter_ingredients(mpc,P)

if ~isempty(mpc.ter_ingredients)
    mpc.recompute_cost_hess = 1;

    mpc.P(:,:) = P;
    mpc.P2(:,:) = 2*P;
end

end
