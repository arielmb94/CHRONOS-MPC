% INIT_MPC_STATE_CNSTR Defines state constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_STATE_CNSTR(mpc, x_min, x_max) sets strict (hard) lower and 
%   upper bounds on the state variables. The solver will strictly enforce these 
%   limits. This is suitable for absolute physical boundaries, but may cause 
%   the solver to crash (go infeasible) if a disturbance pushes the system too far.
%
%   mpc = INIT_MPC_STATE_CNSTR(mpc, x_min, x_max, x_min_slack_active, x_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. Soft constraints 
%   can be safely violated during massive disturbances to keep the solver running, 
%   while applying a customizable penalty to drive the state back within limits 
%   as quickly as possible.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure.
%       x_min              - [nx x 1] Array of lower state limits (use [] if none).
%       x_max              - [nx x 1] Array of upper state limits (use [] if none).
%       x_min_slack_active - (Optional) [nx x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding x_min limit soft.
%       x_max_slack_active - (Optional) [nx x 1] Logical array or scalar. Set to 1 
%                            to make the corresponding x_max limit soft.
%       qv_min             - (Optional) [nx x 1] or scalar. Penalty weight for violating 
%                            the x_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [nx x 1] or scalar. Penalty weight for violating 
%                            the x_max soft limits. Higher values mean stricter enforcement.
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
%         that setting across all constrained states.
%       - Example: Passing x_max_slack_active = [0; 1; 0] makes only the second 
%         state limit soft, leaving the first and third as hard constraints.
function mpc = init_mpc_state_cnstr(mpc,x_min,x_max,...
                x_min_slack_active,x_max_slack_active,qv_min,qv_max)
arguments
    mpc
    x_min
    x_max
    x_min_slack_active = []
    x_max_slack_active = []
    qv_min = []
    qv_max = []
end

s_cnstr.min = x_min;
s_cnstr.max = x_max;

is_there_slack = 0;

if ~isempty(s_cnstr.min)
    
    s_cnstr.fi_min_x0 = zeros(mpc.Nx,1);

    s_cnstr.grad_min = -1 * genGradX(mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(x_min_slack_active)

        is_there_slack = 1;

        [mpc,s_cnstr] = init_slack_min_condition(mpc,s_cnstr,x_min_slack_active,...
                        qv_min,mpc.Nx,mpc.nx);
    else
        s_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [s_cnstr.hess_min,mi] = genHessIneq(s_cnstr.grad_min);
    mpc.m = mpc.m+mi;

end

if ~isempty(s_cnstr.max)
    
    s_cnstr.fi_max_x0 = zeros(mpc.Nx,1);

    s_cnstr.grad_max = genGradX(mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(x_max_slack_active)

        is_there_slack = 1;

        [mpc,s_cnstr] = init_slack_max_condition(mpc,s_cnstr,x_max_slack_active,...
                        qv_max,mpc.Nx,mpc.nx);
    else
        s_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [s_cnstr.hess_max,mi] = genHessIneq(s_cnstr.grad_max);
    mpc.m = mpc.m+mi;
    
end

mpc.s_cnstr = s_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end