% INIT_MPC_OUTPUT_CNSTR Defines output constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_OUTPUT_CNSTR(mpc, y_min, y_max) sets strict (hard) lower and 
%   upper bounds on the output. The solver will strictly enforce these 
%   limits. This is suitable for absolute physical boundaries, but may cause 
%   the solver to crash (go infeasible) if a disturbance pushes the system too far.
%
%   mpc = INIT_MPC_OUTPUT_CNSTR(mpc, y_min, y_max, y_min_slack_active, y_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. Soft constraints 
%   can be safely violated during massive disturbances to keep the solver running, 
%   while applying a customizable penalty to drive the state back within limits 
%   as quickly as possible.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       y_min              - [ny x 1] Array of lower output limits (use [] if none).
%       y_max              - [ny x 1] Array of upper output limits (use [] if none).
%       y_min_slack_active - (Optional) [ny x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding y_min limit soft.
%       y_max_slack_active - (Optional) [ny x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding y_max limit soft.
%       qv_min             - (Optional) [ny x 1] or scalar. Penalty weight for violating 
%                            the y_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [ny x 1] or scalar. Penalty weight for violating 
%                            the y_max soft limits. Higher values mean stricter enforcement.
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
%         that setting across all constrained outputs.
%       - Example: Passing y_max_slack_active = [0; 1; 0] makes only the second 
%         output limit soft, leaving the first and third as hard constraints.
function mpc = init_mpc_output_cnstr(mpc,y_min,y_max,...
                y_min_slack_active,y_max_slack_active,qv_min,qv_max)
arguments
    mpc
    y_min
    y_max
    y_min_slack_active = []
    y_max_slack_active = []
    qv_min = []
    qv_max = []
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(y_min, mpc.ny, 'y_min');
validate_column_vector(y_max, mpc.ny, 'y_max');
validate_column_vector(y_min_slack_active, mpc.ny, 'y_min_slack_active');
validate_column_vector(y_max_slack_active, mpc.ny, 'y_max_slack_active');
validate_column_vector(qv_min, mpc.ny, 'qv_min');
validate_column_vector(qv_max, mpc.ny, 'qv_max');

% Expand scalars to full local vectors if needed
if isscalar(y_min), y_min = y_min * ones(mpc.ny, 1); end
if isscalar(y_max), y_max = y_max * ones(mpc.ny, 1); end

y_cnstr.min = y_min;
y_cnstr.max = y_max;

is_there_slack = 0;

if ~isempty(y_cnstr.min)

    y_cnstr.fi_min_x0 = zeros(mpc.Ny,1);

    % Outputs box constraints
    y_cnstr.grad_min = -1 * genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
                            mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(y_min_slack_active)

        is_there_slack = 1;

        [mpc,y_cnstr] = init_slack_min_condition(mpc,y_cnstr,...
                        y_min_slack_active,qv_min,mpc.Ny,mpc.ny);
    else
        y_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [y_cnstr.hess_min,mi] = genHessIneq(y_cnstr.grad_min);
    mpc.m = mpc.m+mi;
    
end

if ~isempty(y_cnstr.max)

    y_cnstr.fi_max_x0 = zeros(mpc.Ny,1);

    % Outputs box constraints
    y_cnstr.grad_max = genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
                       mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(y_max_slack_active)

        is_there_slack = 1;

        [mpc,y_cnstr] = init_slack_max_condition(mpc,y_cnstr,...
                        y_max_slack_active,qv_max,mpc.Ny,mpc.ny);
    else
        y_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [y_cnstr.hess_max,mi] = genHessIneq(y_cnstr.grad_max);
    mpc.m = mpc.m+mi;

end

mpc.y_cnstr = y_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end