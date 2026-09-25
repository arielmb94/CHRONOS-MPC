% UPDATE_MPC_TRACKING_COST Update the tracking-error weight.
%
%   mpc = UPDATE_MPC_TRACKING_COST(mpc, Qe) replaces the weight in
%
%       J_tracking += 0.5*(r_k-y_k)'*Qe_k*(r_k-y_k).
%
%   Qe may be ny-by-ny or ny-by-ny-by-L, where L is the number of supplied
%   horizon stages. If L < N, the last supplied stage is reused for the
%   remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Qe      - Updated tracking weight, size ny-by-ny or ny-by-ny-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
function mpc = update_mpc_Tracking_cost(mpc,Qe)

if ~isempty(mpc.tracking_cost)
    mpc.update_tracking = true;

    len_Q = size(Qe,3);
    mpc.Qe(:,:,:) = fill_mat(mpc.Qe, Qe, 1);
    if ~isempty(mpc.y_use_ter)
        ter_stage = len_Q;
        if ter_stage > mpc.N, ter_stage = mpc.N; end
        mpc.Qe_ter(:,:) = Qe(mpc.y_rows_ter,mpc.y_rows_ter,ter_stage);
    end

    if ~isempty(mpc.y_use_k0), mpc.Qe_0(:,:) = Qe(mpc.y_rows_k0,mpc.y_rows_k0,1); end
end

end
