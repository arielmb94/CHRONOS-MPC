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
function mpc = update_mpc_ControlRate_cost(mpc,Rdu)

mpc.recompute_cost_hess = 1;

len_Rdu = size(Rdu,3);
if len_Rdu < mpc.N
    mpc.Rdu(:,:,:) = fill_mat(mpc.Rdu, Rdu, 1);
else
    mpc.Rdu(:,:,:) = Rdu(:,:,1:mpc.N);
end

for k = 1:mpc.N-1
    ku = k+1;
    mpc.gradRateCtrl_Rdu_k(:,:,k) = [-mpc.Rdu(:,:,ku);mpc.Rdu(:,:,ku)];
    mpc.H_RateCtrl_k(:,:,k) = [mpc.Rdu(:,:,ku) -mpc.Rdu(:,:,ku);
                               -mpc.Rdu(:,:,ku) mpc.Rdu(:,:,ku)];
end

end