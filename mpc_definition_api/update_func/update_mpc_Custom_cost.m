% UPDATE_MPC_CUSTOM_COST Update custom-cost weights.
%
%   mpc = UPDATE_MPC_CUSTOM_COST(mpc, Qz, qz) replaces the quadratic
%   and/or linear weights applied to the custom signal z_k:
%
%       J_custom += 0.5*z_k'*Qz_k*z_k + qz_k'*z_k.
%
%   Use [] to leave either weight unchanged. The corresponding quadratic
%   or linear term must first be enabled with INIT_MPC_CUSTOM_COST. To
%   update the definition of z_k, use UPDATE_MPC_CUSTOM_COST_VECTOR.
%
%   Qz may be nz-by-nz or nz-by-nz-by-L, and qz may be nz-by-1 or nz-by-L.
%   L is the number of supplied horizon stages. If L < N, the last supplied
%   stage is reused for the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Qz      - Optional updated quadratic weight, size nz-by-nz or
%               nz-by-nz-by-L.
%     qz      - Optional updated linear weight, size nz-by-1 or nz-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the quadratic weight:
%
%       mpc = update_mpc_Custom_cost(mpc, Qz, []);
function mpc = update_mpc_Custom_cost(mpc,Qz,qz)

if ~isempty(Qz) && ~isempty(mpc.quad_custom_cost)
    mpc.update_customcost_quad = true;

    len_Q = size(Qz,3);
    mpc.Qz(:,:,:) = fill_mat(mpc.Qz, Qz, 1);
    if ~isempty(mpc.z_use_ter)
        ter_stage = len_Q;
        if ter_stage > mpc.N, ter_stage = mpc.N; end
        mpc.Qz_ter(:,:) = Qz(mpc.z_rows_ter,mpc.z_rows_ter,ter_stage);
    end

    if ~isempty(mpc.z_use_k0), mpc.Qz_0(:,:) = Qz(mpc.z_rows_k0,mpc.z_rows_k0,1); end
end

if ~isempty(qz) && ~isempty(mpc.lin_custom_cost)
    mpc.update_customcost_lin = true;

    len_q = size(qz,2);
    mpc.qz(:,:) = fill_vec(mpc.qz, qz, 1);
    if ~isempty(mpc.z_use_ter)
        ter_stage = len_q;
        if ter_stage > mpc.N, ter_stage = mpc.N; end
        mpc.qz_ter(:) = qz(mpc.z_rows_ter,ter_stage);
    end

    if ~isempty(mpc.z_use_k0), mpc.qz_0(:) = qz(mpc.z_rows_k0,1); end
end
    
end
