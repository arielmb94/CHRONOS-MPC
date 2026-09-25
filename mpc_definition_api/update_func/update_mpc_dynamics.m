% UPDATE_MPC_DYNAMICS Update the prediction-model matrices.
%
%   mpc = UPDATE_MPC_DYNAMICS(mpc, A, B, Bd) updates selected matrices in
%
%       s_(k+1) = A_k*s_k + B_k*u_k + Bd_k*d_k.
%
%   Use [] to leave a matrix unchanged. The dynamics must first be defined
%   with INIT_MPC_DYNAMICS. Online updates cannot change the matrix dimensions
%   nor add new terms not enabled during initialization.
%
%   A, B, and Bd may be constant matrices or contain L horizon stages in
%   their third dimension. If L < N, the last supplied stage is reused for
%   the remaining stages.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     A       - Optional state matrix, size nx-by-nx or nx-by-nx-by-L.
%     B       - Optional control matrix, size nx-by-nu or nx-by-nu-by-L.
%     Bd      - Optional fixed-known-input matrix, size nx-by-nd or
%               nx-by-nd-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only A and B:
%
%       mpc = update_mpc_dynamics(mpc, A, B, []);
function mpc = update_mpc_dynamics(mpc,A,B,Bd)

if ~isempty(A)
    mpc.A(:,:,:) = fill_mat(mpc.A, A, 1);
end

if ~isempty(B)
    mpc.B(:,:,:) = fill_mat(mpc.B, B, 1);
end

if ~isempty(Bd) && ~isempty(mpc.dyn_use_d)
    mpc.Bd(:,:,:) = fill_mat(mpc.Bd, Bd, 1);
end

end
