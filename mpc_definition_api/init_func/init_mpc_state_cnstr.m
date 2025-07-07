%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_state_cnstr(mpc,x_min,x_max)
%
% Define constraints on the state vector x:
%
%   x_min <= x <= x_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - x_min (optional): nx column vector, lower bound constraint values on 
%   the state vector
%   - x_max (optional): nx column vector, upper bound constraint values on 
%   the state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_state_cnstr(mpc,x_min,x_max)

s_cnstr.min = x_min;
s_cnstr.max = x_max;

[s_cnstr.grad_min,s_cnstr.grad_max] = genGradX(mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(s_cnstr.min)
    s_cnstr.fi_min_x0 = zeros(mpc.Nx,1);
    [s_cnstr.hess_min,mi] = genHessIneq(s_cnstr.grad_min);
    mpc.m = mpc.m+mi;
end
if ~isempty(s_cnstr.max)
    s_cnstr.fi_max_x0 = zeros(mpc.Nx,1);
    [s_cnstr.hess_max,mi] = genHessIneq(s_cnstr.grad_max);
    mpc.m = mpc.m+mi;
end

mpc.s_cnstr = s_cnstr;
end