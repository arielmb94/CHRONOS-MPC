%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_Tracking_cost(mpc,Qe)
%
% Adds quadratic penalties on the tracking error:
%
%   J += (r - y)' * Qe * (r - y)
%
% y is the tracking feedback signal, defined during the call to 
% init_mpc_system(), the reference vector r is introduced during MPC 
% runtime iterations on the call to mpc_solve().
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Qe: ny x ny square matrix, weights for the quadratic penalty on the
%   tracking error
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_Tracking_cost(mpc,Qe)

mpc.Qe = Qe;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
end

[mpc.gradErrQe,mpc.hessErrTerm] = genLinOutGradHess(Qe,mpc.C,mpc.D,mpc.N,...
        mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

mpc.hessCost = mpc.hessCost + mpc.hessErrTerm;

end