%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    mpc = init_mpc_DiffControl_cost(mpc,Rdu)
%
% Adds penalty on the control variation between sampling instances:
%
%   J += delta_u' * Rdu * delta_u
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
function mpc = init_mpc_DiffControl_cost(mpc,Rdu)

mpc.Rdu = Rdu;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx+mpc.Nv);
end

mpc.gradDiffCtlrR = zeros(mpc.Nu+mpc.Nx+mpc.Nv,mpc.Nu);
mpc.gradDiffCtlr = zeros(mpc.Nu+mpc.Nx+mpc.Nv,mpc.Nu);
mpc.hessDiffCtrlTerm = zeros(mpc.Nu+mpc.Nx+mpc.Nv);

mpc = genDiffControlGradHess(mpc,Rdu,mpc.N_ctr_hor,mpc.nx,mpc.nu);

mpc.hessCost = mpc.hessCost + mpc.hessDiffCtrlTerm;

end