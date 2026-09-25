% UPDATE_MPC_CUSTOM_CNSTR_VECTOR Update the custom-constraint signal model.
%
%   mpc = UPDATE_MPC_CUSTOM_CNSTR_VECTOR(mpc, Ch, Dh, Dsuh, Ddh) updates
%   selected coefficients in
%
%       h_k = Ch_k*s_k + Dh_k*u_k + Dsuh_k*su_k + Ddh_k*dh_k.
%
%   Use [] to leave a coefficient unchanged. The signal and its bounds must
%   first be defined with INIT_MPC_CUSTOM_CNSTR. This updater cannot change
%   the number of custom constraints nor add a coefficient term that was not
%   enabled during initialization.
%
%   Coefficients may be constant matrices or contain L horizon stages in
%   their third dimension. If L < N, the last supplied stage is reused for
%   the remaining stages. Here, su_k is the control action preceding u_k,
%   and dh_k is the dedicated fixed known input supplied to MPC_SOLVE.
% 
%   Existing custom bounds use the updated signal model. To
%   change the bounds, use UPDATE_MPC_CUSTOM_CNSTR_LIMITS.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     Ch      - Optional state coefficient, size nh-by-nx or nh-by-nx-by-L.
%     Dh      - Optional control coefficient, size nh-by-nu or
%               nh-by-nu-by-L.
%     Dsuh    - Optional previous-control coefficient, size nh-by-nu or
%               nh-by-nu-by-L.
%     Ddh     - Optional fixed-known-input coefficient, size nh-by-ndh or
%               nh-by-ndh-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only Ch and Dh:
%
%       mpc = update_mpc_custom_cnstr_vector(mpc, Ch, Dh, [], []);
function mpc = update_mpc_custom_cnstr_vector(mpc,Ch,Dh,Dsuh,Ddh)

if ~isempty(Ch) && ~isempty(mpc.h_cnstr.use_s)

    len_Ch = size(Ch,3);
    mpc.Ch(:,:,:) = fill_mat(mpc.Ch, Ch, 1);
    if ~isempty(mpc.h_cnstr.use_ter)
        ter_stage = len_Ch;
        if ter_stage > mpc.N, ter_stage = mpc.N; end
        mpc.Ch_ter(:,:) = Ch(mpc.h_cnstr.rows_ter,:,ter_stage);
    end
    if ~isempty(mpc.h_cnstr.use_k0), mpc.Ch_0(:,:) = Ch(mpc.h_cnstr.rows_k0,:,1); end

end

if ~isempty(Dh) && ~isempty(mpc.h_cnstr.use_u)
    mpc.Dh(:,:,:) = fill_mat(mpc.Dh, Dh, 1);
    if ~isempty(mpc.h_cnstr.use_k0)
        mpc.Dh_0(:,:) = Dh(mpc.h_cnstr.rows_k0,:,1);
    end

end

if ~isempty(Dsuh) && ~isempty(mpc.h_cnstr.use_su)
    mpc.Dsuh(:,:,:) = fill_mat(mpc.Dsuh, Dsuh, 1);
    if ~isempty(mpc.h_cnstr.use_k0), mpc.Dsuh_0(:,:) = Dsuh(mpc.h_cnstr.rows_k0,:,1); end

end

if ~isempty(Ddh) && ~isempty(mpc.h_cnstr.use_d)
    mpc.Ddh(:,:,:) = fill_mat(mpc.Ddh, Ddh, 1);
    if ~isempty(mpc.h_cnstr.use_k0), mpc.Ddh_0(:,:) = Ddh(mpc.h_cnstr.rows_k0,:,1); end
end

end
