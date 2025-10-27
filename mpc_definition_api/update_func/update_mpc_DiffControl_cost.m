%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_DiffControl_cost(mpc,Rdu)
%
% Modifies the weight Rdu for the quadratic penalty term on the control 
% action variation between sampling instances. The function then updates 
% the MPC gradients and Hessians accordingly.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Rdu: nu x nu square matrix, weights for the quadratic penalty term on
%   the control control variation between sampling instances
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_DiffControl_cost(mpc,Rdu)

mpc.Rdu = Rdu;

mpc = genDiffControlGradHess(mpc,Rdu,mpc.N_ctr_hor,mpc.nx,mpc.nu);

mpc.recompute_cost_hess = 1;

end