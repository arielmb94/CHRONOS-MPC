% INIT_MPC_U_CNSTR Defines control action constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_U_CNSTR(mpc, u_min, u_max) sets strict (hard) lower and 
%   upper bounds on the control action. The solver will strictly enforce these 
%   limits.
%
%   mpc = INIT_MPC_U_CNSTR(mpc, u_min, u_max, u_min_slack_active, u_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. For control
%   actions it is advisable to use hard constraints since it is a magnitude
%   the solver will have direct control over.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       u_min              - [nu x 1] Array of lower control action limits (use [] if none).
%       u_max              - [nu x 1] Array of upper control action limits (use [] if none).
%       u_min_slack_active - (Optional) [nu x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding u_min limit soft.
%       u_max_slack_active - (Optional) [nu x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding u_max limit soft.
%       qv_min             - (Optional) [nu x 1] or scalar. Penalty weight for violating 
%                            the u_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [nu x 1] or scalar. Penalty weight for violating 
%                            the u_max soft limits. Higher values mean stricter enforcement.
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (constraint gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
%
%   USAGE TIPS:
%       - If soft constraitns are enabled and qv_min or qv_max are not
%         passed, the soft constraint penalty weight will default to the value 
%         stored in mpc.qv
%       - Passing a scalar to the slack or qv inputs will automatically apply 
%         that setting across all constrained control actions.
%       - Example: Passing u_max_slack_active = [0; 1; 0] makes only the second 
%         control action limit soft, leaving the first and third as hard constraints.
function mpc = init_mpc_u_cnstr(mpc,u_min,u_max,u_min_slack_active,u_max_slack_active,...
                                qv_min,qv_max)
arguments
    mpc
    u_min
    u_max
    u_min_slack_active = []
    u_max_slack_active = []
    qv_min = []
    qv_max = []
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(u_min, mpc.nu, 'u_min');
validate_column_vector(u_max, mpc.nu, 'u_max');
validate_column_vector(u_min_slack_active, mpc.nu, 'u_min_slack_active');
validate_column_vector(u_max_slack_active, mpc.nu, 'u_max_slack_active');
validate_column_vector(qv_min, mpc.nu, 'qv_min');
validate_column_vector(qv_max, mpc.nu, 'qv_max');

% Expand scalars to full local vectors if needed
if isscalar(u_min), u_min = u_min * ones(mpc.nu, 1); end
if isscalar(u_max), u_max = u_max * ones(mpc.nu, 1); end

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
                        qv_min,mpc.Nu,mpc.nu);
    else
        u_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [u_cnstr.hess_min,mi] = genHessIneq(u_cnstr.grad_min);
    mpc.m = mpc.m+mi;

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
                        qv_max,mpc.Nu,mpc.nu);
    else
        u_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [u_cnstr.hess_max,mi] = genHessIneq(u_cnstr.grad_max);
    mpc.m = mpc.m+mi;
    
end

mpc.u_cnstr = u_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end