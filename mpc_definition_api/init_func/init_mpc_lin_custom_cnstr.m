% INIT_MPC_LIN_CUSTOM_CNSTR Defines constraints and soft-penalties on
% custom user defined signals.
%
%   mpc = INIT_MPC_LIN_CUSTOM_CNSTR(mpc, Ch, Dh, Ddh, h_min, h_max) defines a 
%   custom auxiliary signal 'h' and sets strict (hard) lower and upper bounds 
%   on it. The custom signal is calculated as:
%
%       h = Ch * x + Dh * u + Ddh * dh
%
%   The solver will strictly enforce h_min <= h <= h_max. This is suitable for 
%   hard physical limits, but massive disturbances could cause the solver to crash 
%   (go infeasible) if the limits are impossible to satisfy.
%
%   mpc = INIT_MPC_LIN_CUSTOM_CNSTR(mpc, ..., h_min_slack_active, h_max_slack_active, qv_min, qv_max) 
%   allows you to define these custom bounds as "soft" constraints. Soft constraints 
%   can be safely violated during severe disturbances to prevent solver crashes, 
%   while applying a customizable linear penalty to drive the signal back within 
%   limits as quickly as actuator power allows. It is highly advisible to
%   enable soft constraint on the custom signals.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure.
%       Ch                 - [nh x nx] Matrix mapping states to the custom signal.
%       Dh                 - [nh x nu] Matrix mapping inputs to the custom signal.
%       Ddh                - [nh x ndh] Matrix mapping measured disturbances to the custom signal.
%       h_min              - [nh x 1] Array of lower limits (use [] if none).
%       h_max              - [nh x 1] Array of upper limits (use [] if none).
%       h_min_slack_active - (Optional) [nh x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding h_min limit soft.
%       h_max_slack_active - (Optional) [nh x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding h_max limit soft.
%       qv_min             - (Optional) [nh x 1] or scalar. Penalty weight for violating 
%                            the h_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [nh x 1] or scalar. Penalty weight for violating 
%                            the h_max soft limits. Higher values mean stricter enforcement.
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (matrices, gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
%
%   EXAMPLE USE CASE: Tracking an input reference (u_star)
%       We want to limit the variation of the control action with respect to a 
%       dynamic target value, meaning we want to constrain: h = u - u_star.
%       To achieve this, we set up our custom signal matrices as:
%           - Ch  = zeros(nu, nx)
%           - Dh  = eye(nu)
%           - Ddh = -eye(nu)
%       This creates the equation: h = 0*x + 1*u - 1*dh.
%       During runtime, the user passes 'u_star' into the 'dh' disturbance 
%       vector, and the solver handles the rest!
%
%   USAGE TIPS:
%       - If soft constraitns are enabled and qv_min or qv_max are not
%         passed, the soft constraint penalty weight will default to the value 
%         stored in mpc.qv
%       - Passing a scalar to the slack or qv inputs will automatically apply 
%         that setting across all constrained outputs.
%       - Example: Passing h_max_slack_active = [0; 1; 0] makes only the second 
%         output limit soft, leaving the first and third as hard constraints.
function mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max, ...
                h_min_slack_active,h_max_slack_active,qv_min,qv_max)
arguments
    mpc
    Ch
    Dh
    Ddh
    h_min
    h_max
    h_min_slack_active = []
    h_max_slack_active = []
    qv_min = []
    qv_max = []
end

% General Inequality Matrix
mpc.Ch = Ch;
mpc.Dh = Dh;
mpc.Ddh = Ddh;

mpc.ndh = size(Ddh,2);  %number of disturbance inputs to general inequalities
mpc.nh = size(Ch,1);  %number of general inequalities

% INPUT DIMENSION VALIDATION 
validate_column_vector(h_min, mpc.nh, 'h_min');
validate_column_vector(h_max, mpc.nh, 'h_max');
validate_column_vector(h_min_slack_active, mpc.nh, 'h_min_slack_active');
validate_column_vector(h_max_slack_active, mpc.nh, 'h_max_slack_active');
validate_column_vector(qv_min, mpc.nh, 'qv_min');
validate_column_vector(qv_max, mpc.nh, 'qv_max');

% Expand scalars to full local vectors if needed
if isscalar(h_min), h_min = h_min * ones(mpc.nh, 1); end
if isscalar(h_max), h_max = h_max * ones(mpc.nh, 1); end

if mpc.Dh == 0
    mpc.Nh = (mpc.N-1)*mpc.nh;
else
    mpc.Nh = mpc.N*mpc.nh;
end

mpc.h = zeros(mpc.Nh,1);

mpc.Ndh = mpc.N*mpc.ndh;

% General Inequalites box constraints
h_cnstr.min = h_min;
h_cnstr.max = h_max;

is_there_slack = 0;

if ~isempty(h_cnstr.min)

    h_cnstr.fi_min_x0 = zeros(mpc.Nh,1);

    h_cnstr.grad_min = -1 * genGradY(mpc.Ch,mpc.Dh,mpc.N,mpc.N_ctr_hor,...
                            mpc.Nx,mpc.Nu,mpc.Nh,mpc.nx,mpc.nu,mpc.nh,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(h_min_slack_active)

        is_there_slack = 1;

        [mpc,h_cnstr] = init_slack_min_condition(mpc,h_cnstr,...
                        h_min_slack_active,qv_min,mpc.Nh,mpc.nh);
    else
        h_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [h_cnstr.hess_min,mi] = genHessIneq(h_cnstr.grad_min);
    mpc.m = mpc.m+mi;

end
if ~isempty(h_cnstr.max)

    h_cnstr.fi_max_x0 = zeros(mpc.Nh,1);

    h_cnstr.grad_max = genGradY(mpc.Ch,mpc.Dh,mpc.N,mpc.N_ctr_hor,...
                       mpc.Nx,mpc.Nu,mpc.Nh,mpc.nx,mpc.nu,mpc.nh,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(h_max_slack_active)

        is_there_slack = 1;

        [mpc,h_cnstr] = init_slack_max_condition(mpc,h_cnstr,...
                        h_max_slack_active,qv_max,mpc.Nh,mpc.nh);
    else
        h_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [h_cnstr.hess_max,mi] = genHessIneq(h_cnstr.grad_max);
    mpc.m = mpc.m+mi;

end

mpc.h_cnstr = h_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end