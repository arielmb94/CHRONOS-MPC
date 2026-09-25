% UPDATE_MPC_OUTPUT_VECTOR Update the tracking-output model.
%
%   mpc = UPDATE_MPC_OUTPUT_VECTOR(mpc, C, D, Dd) updates selected
%   coefficients in
%
%       y_k = C_k*s_k + D_k*u_k + Dd_k*d_k.
%
%   Use [] to leave a coefficient unchanged. The output must first be
%   defined with INIT_MPC_DYNAMICS or INIT_MPC_OUTPUT. This updater cannot
%   change the number of outputs or add a coefficient term that was not
%   enabled during initialization.
%
%   C, D, and Dd may be constant matrices or contain L horizon stages in
%   their third dimension. If L < N, the last supplied stage is reused for
%   the remaining stages. The input d_k is the same fixed known input used
%   by the dynamics.
% 
%   Existing tracking costs and output constraints use the
%   updated output model.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     C       - Optional state coefficient, size ny-by-nx or ny-by-nx-by-L.
%     D       - Optional control coefficient, size ny-by-nu or
%               ny-by-nu-by-L.
%     Dd      - Optional fixed-known-input coefficient, size ny-by-nd or
%               ny-by-nd-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the state and control coefficients:
%
%       mpc = update_mpc_output_vector(mpc, C, D, []);
function mpc = update_mpc_output_vector(mpc,C,D,Dd)

if ~isempty(C) && ~isempty(mpc.y_use_s)
    if ~isempty(mpc.tracking_cost)
        mpc.update_tracking = true;
    end

    len_C = size(C,3);
    mpc.C(:,:,:) = fill_mat(mpc.C, C, 1);
    if ~isempty(mpc.y_use_ter)
        ter_stage = len_C;
        if ter_stage > mpc.N, ter_stage = mpc.N; end
        mpc.C_ter(:,:) = C(mpc.y_rows_ter,:,ter_stage);
    end
    if ~isempty(mpc.y_use_k0), mpc.C_0(:,:) = C(mpc.y_rows_k0,:,1); end

end

if ~isempty(D) && ~isempty(mpc.y_use_u)
    if ~isempty(mpc.tracking_cost)
        mpc.update_tracking = true;
    end
    mpc.D(:,:,:) = fill_mat(mpc.D, D, 1);
    if ~isempty(mpc.y_use_k0)
        mpc.D_0(:,:) = D(mpc.y_rows_k0,:,1);
    end
end

if ~isempty(Dd) && ~isempty(mpc.y_use_d)
    mpc.Dd(:,:,:) = fill_mat(mpc.Dd, Dd, 1);

    if ~isempty(mpc.y_use_k0), mpc.Dd_0(:,:) = Dd(mpc.y_rows_k0,:,1); end
end

end
