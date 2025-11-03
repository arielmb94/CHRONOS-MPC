%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_u_cnstr(mpc,u_min,u_max)
%
% Define limits on the control action:
%
%   u_min <= u <= u_max
%
% You can define either lower or upper bounds alone, or both.
%
% Example uses:
%
%   - only lower bound: mpc = init_mpc_u_cnstr(mpc,u_min,[])
%   - only upper bound: mpc = init_mpc_u_cnstr(mpc,[],u_max)
%   - upper and lower bounds: mpc = init_mpc_u_cnstr(mpc,u_min,u_max)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_u_cnstr(mpc,u_min,u_max,u_min_activ,u_max_activ)

u_cnstr.min = u_min;
u_cnstr.max = u_max;

% Control box constraints
[u_cnstr.grad_min,u_cnstr.grad_max] = genGradU(mpc.N_ctr_hor,...
                                    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(u_cnstr.min)
    u_cnstr.fi_min_x0 = zeros(mpc.Nu,1);
    [u_cnstr.hess_min,mi] = genHessIneq(u_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    u_cnstr.grad_min_feas_slv = [u_cnstr.grad_min;-ones(1,mpc.Nu)];
    
    u_cnstr.hess_min_feas_slv = genHessIneq(u_cnstr.grad_min_feas_slv);

    % Active-set like optimization setting
    if exist("u_min_activ")
        u_cnstr = set_active_set_min(u_cnstr,u_min,u_min_activ,mpc.nu);
    else
        u_cnstr = empty_active_set_min(u_cnstr);
    end
end

if ~isempty(u_cnstr.max)
    u_cnstr.fi_max_x0 = zeros(mpc.Nu,1);
    [u_cnstr.hess_max,mi] = genHessIneq(u_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    u_cnstr.grad_max_feas_slv = [u_cnstr.grad_max;-ones(1,mpc.Nu)];
    
    u_cnstr.hess_max_feas_slv = genHessIneq(u_cnstr.grad_max_feas_slv);

    % Active-set like optimization setting
    if exist("u_max_activ")
        u_cnstr = set_active_set_max(u_cnstr,u_max,u_max_activ,mpc.nu);
    else
        u_cnstr = empty_active_set_max(u_cnstr);
    end
end

mpc.u_cnstr = u_cnstr;
end