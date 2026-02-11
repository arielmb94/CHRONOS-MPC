%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max,...
%                               du_min_slack_active,du_max_slack_active,...
%                               du_min_hard,du_max_hard)
%
% Define constraint on the control variation between sampling instances:
%
%   du_min <= delta_u <= du_max
%
% In some scenarios, it is ok if these limits are violated by a small
% amount. In these cases we can enable slack variables v_i on the
% constraints to allow for violations of the now soft constraint limits,
% the constraint now becoming:
%
%   du_min - delta_u <= v_i
%   delta_u - du_max <= v_i
% 
% Even if small violations on the soft constraints can be tolerated,
% sometimes there are phisical limits which cannot be surpased. In those 
% cases we can enable hard limits such that:
%
% du_min_hard <= du_min <= delta_u <= du_max <= du_max_hard
%
% In:
%   - mpc: CHRONOS mpc structure
%
%   - du_min (optional): nu column vector, lower bound constraint values on
%   the control variation between sampling instances
%
%   - du_max (optional): nu column vector, upper bound constraint values on
%   the control variation between sampling instances
%
%   - du_min_slack_active (optional): single boolean or nu boolean column 
%   vector, indicates which elements of the differential control vector 
%   minumum constraints have slack variables
%
%   - du_max_slack_active (optional): single boolean or nu boolean column 
%   vector, indicates which elements of the differential control vector 
%   maximum constraints have slack variables
%
%   - du_min_hard (optional): nu column vector, maximum lower bound 
%   constraint values on the differential control vector
%
%   - du_max_hard (optional): nu column vector, maximum upper bound 
%   constraint values on the differential control vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max,...
               du_min_slack_active,du_max_slack_active,...
               du_min_hard,du_max_hard)
arguments
    mpc
    du_min
    du_max
    du_min_slack_active = []
    du_max_slack_active = []
    du_min_hard = []
    du_max_hard = []
end

du_cnstr.min = du_min;
du_cnstr.max = du_max;

is_there_slack = 0;

if ~isempty(du_cnstr.min)

    du_cnstr.fi_min_x0 = zeros(mpc.Nu,1);

    % Differential Control box constraints
    du_cnstr.grad_min = -1 * genGradDeltaU(mpc.N_ctr_hor,...
                             mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(du_min_slack_active)

        is_there_slack = 1;

        [mpc,du_cnstr] = init_slack_min_condition(mpc,du_cnstr,...
                        du_min_slack_active,du_min_hard,mpc.Nu,mpc.nu);
    else
        du_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [du_cnstr.hess_min,mi] = genHessIneq(du_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    du_cnstr.grad_min_feas_slv = [du_cnstr.grad_min;-ones(1,mpc.Nu)];
    
    du_cnstr.hess_min_feas_slv = genHessIneq(du_cnstr.grad_min_feas_slv);

end

if ~isempty(du_cnstr.max)
    
    du_cnstr.fi_max_x0 = zeros(mpc.Nu,1);

    % Differential Control box constraints
    du_cnstr.grad_max = genGradDeltaU(mpc.N_ctr_hor,...
                        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(du_max_slack_active)

        is_there_slack = 1;

        [mpc,du_cnstr] = init_slack_max_condition(mpc,du_cnstr,...
                        du_max_slack_active,du_max_hard,mpc.Nu,mpc.nu);
    else
        du_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [du_cnstr.hess_max,mi] = genHessIneq(du_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    du_cnstr.grad_max_feas_slv = [du_cnstr.grad_max;-ones(1,mpc.Nu)];
    
    du_cnstr.hess_max_feas_slv = genHessIneq(du_cnstr.grad_max_feas_slv);

end

mpc.du_cnstr = du_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end