%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max)
%
% Define constraint on the control variation between sampling instances:
%
%   du_min <= delta_u <= du_max
%
% You can define either lower or upper bounds alone, or both.
%
% Example uses:
%
%   - only lower bound: mpc = init_mpc_delta_u_cnstr(mpc,du_min,[])
%   - only upper bound: mpc = init_mpc_delta_u_cnstr(mpc,[],du_max)
%   - upper and lower bounds: mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - du_min (optional): nu column vector, lower bound constraint values on
%   the control variation between sampling instances
%   - du_max (optional): nu column vector, upper bound constraint values on
%   the control variation between sampling instances
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max,du_min_activ,du_max_activ)

du_cnstr.min = du_min;
du_cnstr.max = du_max;

% Differential Control box constraints
[du_cnstr.grad_min,du_cnstr.grad_max] = genGradDeltaU(mpc.N_ctr_hor,...
                                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(du_cnstr.min)
    du_cnstr.fi_min_x0 = zeros(mpc.Nu,1);
    [du_cnstr.hess_min,mi] = genHessIneq(du_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    du_cnstr.grad_min_feas_slv = [du_cnstr.grad_min;-ones(1,mpc.Nu)];
    
    du_cnstr.hess_min_feas_slv = genHessIneq(du_cnstr.grad_min_feas_slv);

    % Active-set like optimization setting
    if exist("du_min_activ")
        du_cnstr = set_active_set_min(du_cnstr,du_min,du_min_activ,mpc.nu);
    else
        du_cnstr = empty_active_set_min(du_cnstr);
    end

end

if ~isempty(du_cnstr.max)
    du_cnstr.fi_max_x0 = zeros(mpc.Nu,1);
    [du_cnstr.hess_max,mi] = genHessIneq(du_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    du_cnstr.grad_max_feas_slv = [du_cnstr.grad_max;-ones(1,mpc.Nu)];
    
    du_cnstr.hess_max_feas_slv = genHessIneq(du_cnstr.grad_max_feas_slv);

    % Active-set like optimization setting
    if exist("du_max_activ")
        du_cnstr = set_active_set_max(du_cnstr,du_max,du_max_activ,mpc.nu);
    else
        du_cnstr = empty_active_set_max(du_cnstr);
    end

end

mpc.du_cnstr = du_cnstr;

end