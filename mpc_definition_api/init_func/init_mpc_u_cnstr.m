%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_u_cnstr(mpc,u_min,u_max,...
%                          u_min_slack_active,u_max_slack_active,...
%                          u_min_hard,u_max_hard)
%
% Define limits on the control action:
%
%   u_min <= u <= u_max
%
% In some scenarios, it is ok if these limits are violated by a small
% amount. In these cases we can enable slack variables v_i on the
% constraints to allow for violations of the now soft constraint limits,
% the constraint now becoming:
%
%   u_min - u <= v_i
%   u - u_max <= v_i
% 
% Even if small violations on the soft constraints can be tolerated,
% sometimes there are phisical limits which cannot be surpased. In those 
% cases we can enable hard limits such that:
%
% u_min_hard <= u_min <= u <= u_max <= u_max_hard
%
% In:
%   - mpc: CHRONOS mpc structure
%
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action
%
%   - u_min_slack_active (optional): single boolean or nu boolean column 
%   vector, indicates which elements of the control vector minumum constraints
%   have slack variables
%
%   - u_max_slack_active (optional): single boolean or nu boolean column 
%   vector, indicates which elements of the control vector maximum constraints
%   have slack variables
%
%   - u_min_hard (optional): nu column vector, maximum lower bound 
%   constraint values on the control vector
%
%   - u_max_hard (optional): nu column vector, maximum upper bound 
%   constraint values on the control vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_u_cnstr(mpc,u_min,u_max,...
                u_min_slack_active,u_max_slack_active,u_min_hard,u_max_hard)
arguments
    mpc
    u_min
    u_max
    u_min_slack_active = []
    u_max_slack_active = []
    u_min_hard = []
    u_max_hard = []
end

u_cnstr.min = u_min;
u_cnstr.max = u_max;

is_there_slack = 0;

if ~isempty(u_cnstr.min)
    u_cnstr.fi_min_x0 = zeros(mpc.Nu,1);

    % Control box constraints
    u_cnstr.grad_min = -1 * genGradU(mpc.N_ctr_hor,...
                                     mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(u_min_slack_active)

        is_there_slack = 1;

        [mpc,u_cnstr] = init_slack_min_condition(mpc,u_cnstr,u_min_slack_active,...
                        u_min_hard,mpc.Nu,mpc.nu);
    else
        u_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [u_cnstr.hess_min,mi] = genHessIneq(u_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    u_cnstr.grad_min_feas_slv = [u_cnstr.grad_min;-ones(1,mpc.Nu)];
    
    u_cnstr.hess_min_feas_slv = genHessIneq(u_cnstr.grad_min_feas_slv);

end

if ~isempty(u_cnstr.max)
    u_cnstr.fi_max_x0 = zeros(mpc.Nu,1);

    % Control box constraints
    u_cnstr.grad_max = genGradU(mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(u_max_slack_active)

        is_there_slack = 1;

        [mpc,u_cnstr] = init_slack_max_condition(mpc,u_cnstr,u_max_slack_active,...
                        u_max_hard,mpc.Nu,mpc.nu);
    else
        u_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [u_cnstr.hess_max,mi] = genHessIneq(u_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    u_cnstr.grad_max_feas_slv = [u_cnstr.grad_max;-ones(1,mpc.Nu)];
    
    u_cnstr.hess_max_feas_slv = genHessIneq(u_cnstr.grad_max_feas_slv);

end

mpc.u_cnstr = u_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end