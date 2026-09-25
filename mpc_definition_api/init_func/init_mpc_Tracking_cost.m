% INIT_MPC_TRACKING_COST Add a quadratic tracking-error penalty.
%
%   mpc = INIT_MPC_TRACKING_COST(mpc, Qe) penalizes the error between the
%   runtime reference and the configured tracking output:
%
%       err_k = r_k - y_k
%       J_tracking += 0.5 * err_k' * Qe_k * err_k
%
%   Define the model first with INIT_MPC_DYNAMICS. By default, its tracking
%   output is the full state; call INIT_MPC_OUTPUT before this function to
%   use a different output. Supply the reference as r_in when calling
%   MPC_SOLVE.
%
%   Qe may contain one matrix or L horizon stages, where L is the number of
%   supplied stages. If L < mpc.N, CHRONOS reuses the last stage for the
%   remaining stages, including the terminal stage. If L >= mpc.N, stage
%   mpc.N is used at the terminal stage and later stages are ignored.
%
%   Call this function after defining the dynamics and optional tracking
%   output, and before calling BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     Qe      - Tracking-error weight, size ny-by-ny or ny-by-ny-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - use the same output weight at every stage:
%
%       Qe = diag([10, 1]);
%       mpc = init_mpc_Tracking_cost(mpc, Qe);
function mpc = init_mpc_Tracking_cost(mpc,Qe)

if isempty(Qe) || ~any(Qe(:))
    return;
end

validate_matrix(Qe, mpc.ny, mpc.ny, 'Qe', true);
if isscalar(Qe)
    Qe = Qe * eye(mpc.ny);
elseif size(Qe,1) == 1 && size(Qe,2) == 1
    Qe_staged = Qe;
    Qe = zeros(mpc.ny, mpc.ny, size(Qe_staged,3));
    Qe_eye = eye(mpc.ny);
    for k = 1:size(Qe_staged,3)
        Qe(:,:,k) = Qe_staged(1,1,k) * Qe_eye;
    end
end

mpc.tracking_cost = 1;

mpc.Qe = zeros(mpc.ny,mpc.ny,mpc.N-1);
len_Qe = size(Qe,3);
if len_Qe < mpc.N
    mpc.Qe = fill_mat(mpc.Qe, Qe, 1);
    Qe_ter = mpc.Qe(:,:,mpc.N-1);
else
    mpc.Qe = Qe(:,:,1:mpc.N-1);
    Qe_ter = Qe(:,:,mpc.N);
end

if ~isempty(mpc.y_use_k0)
    mpc.Qe_0 = Qe(mpc.y_rows_k0,mpc.y_rows_k0,1);
else
    mpc.Qe_0 = [];
end
if ~isempty(mpc.y_use_ter)
    mpc.Qe_ter = Qe_ter(mpc.y_rows_ter,mpc.y_rows_ter);
else
    mpc.Qe_ter = [];
end

end
