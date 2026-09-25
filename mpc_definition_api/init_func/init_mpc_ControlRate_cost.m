% INIT_MPC_CONTROLRATE_COST Add a quadratic control-rate penalty.
%
%   mpc = INIT_MPC_CONTROLRATE_COST(mpc, Rdu) penalizes changes in the
%   control action:
%
%       delta_u_0 = u_0 - u_prev
%       delta_u_k = u_k - u_(k-1)
%       J_rate += 0.5 * delta_u_k' * Rdu_k * delta_u_k
%
%   Rdu may contain one matrix or L horizon stages, where L is the number of
%   supplied stages. If L < mpc.N, CHRONOS reuses the last stage for the
%   remaining stages; if L >= mpc.N, only the first mpc.N stages are used.
%
%   Call this function after INIT_MPC_DYNAMICS and before
%   BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     Rdu     - Control-rate weight, size nu-by-nu or nu-by-nu-by-L. Each
%               active stage must be symmetric positive semidefinite.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - use the same rate penalty at every stage:
%
%       Rdu = diag([1, 0.5]);
%       mpc = init_mpc_ControlRate_cost(mpc, Rdu);
function mpc = init_mpc_ControlRate_cost(mpc,Rdu)

if isempty(Rdu) || ~any(Rdu(:))
    return;
end

validate_matrix(Rdu, mpc.nu, mpc.nu, 'Rdu', true);
if isscalar(Rdu)
    Rdu = Rdu * eye(mpc.nu);
elseif size(Rdu,1) == 1 && size(Rdu,2) == 1
    Rdu_staged = Rdu;
    Rdu = zeros(mpc.nu, mpc.nu, size(Rdu_staged,3));
    Rdu_eye = eye(mpc.nu);
    for k = 1:size(Rdu_staged,3)
        Rdu(:,:,k) = Rdu_staged(1,1,k) * Rdu_eye;
    end
end

mpc.has_du = 1;
mpc.controlrate_cost = 1;

mpc.Rdu = zeros(mpc.nu,mpc.nu,mpc.N);
len_Rdu = size(Rdu,3);
if len_Rdu < mpc.N
    mpc.Rdu = fill_mat(mpc.Rdu, Rdu, 1);
else
    mpc.Rdu = Rdu(:,:,1:mpc.N);
end

end
