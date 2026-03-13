% INIT_MPC_DELTA_U_CNSTR Defines control action rate constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_DELTA_U_CNSTR(mpc, du_min, du_max) sets strict (hard) lower and 
%   upper bounds on the control action rate. The solver will strictly enforce these 
%   limits.
%
%   mpc = INIT_MPC_DELTA_U_CNSTR(mpc, du_min, du_max, du_min_slack_active, du_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. For control
%   actions it is advisable to use hard constraints since it is a magnitude
%   the solver will have direct control over.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       du_min              - [nu x 1] Array of lower control action rate limits (use [] if none).
%       du_max              - [nu x 1] Array of upper control action rate limits (use [] if none).
%       du_min_slack_active - (Optional) [nu x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding du_min limit soft.
%       du_max_slack_active - (Optional) [nu x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding du_max limit soft.
%       qv_min             - (Optional) [nu x 1] or scalar. Penalty weight for violating 
%                            the du_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [nu x 1] or scalar. Penalty weight for violating 
%                            the du_max soft limits. Higher values mean stricter enforcement.
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
%         control action rate limit soft, leaving the first and third as hard constraints.
function mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max,...
               du_min_slack_active,du_max_slack_active,qv_min,qv_max)
arguments
    mpc
    du_min
    du_max
    du_min_slack_active = []
    du_max_slack_active = []
    qv_min = []
    qv_max = []
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(du_min, mpc.nu, 'du_min');
validate_column_vector(du_max, mpc.nu, 'du_max');
validate_column_vector(du_min_slack_active, mpc.nu, 'du_min_slack_active');
validate_column_vector(du_max_slack_active, mpc.nu, 'du_max_slack_active');
validate_column_vector(qv_min, mpc.nu, 'qv_min');
validate_column_vector(qv_max, mpc.nu, 'qv_max');

% Expand scalars to full local vectors if needed
if isscalar(du_min), du_min = du_min * ones(mpc.nu, 1); end
if isscalar(du_max), du_max = du_max * ones(mpc.nu, 1); end

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
                        du_min_slack_active,qv_min,mpc.Nu,mpc.nu);
    else
        du_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [du_cnstr.hess_min,mi] = genHessIneq(du_cnstr.grad_min);
    mpc.m = mpc.m+mi;

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
                        du_max_slack_active,qv_max,mpc.Nu,mpc.nu);
    else
        du_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [du_cnstr.hess_max,mi] = genHessIneq(du_cnstr.grad_max);
    mpc.m = mpc.m+mi;

end

mpc.du_cnstr = du_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end