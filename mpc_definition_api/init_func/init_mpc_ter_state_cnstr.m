%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max,...
%                           x_ter_min_slack_active,x_ter_max_slack_active,...
%                           x_ter_min_hard,x_ter_max_hard)
%
% Define box constraints only on the final state vector xN of the MPC 
% prediction horizon:
%
%   xN_min <= xN <= xN_max
%
% In some scenarios, it is ok if these limits are violated by a small
% amount. In these cases we can enable slack variables v_i on the
% constraints to allow for violations of the now soft constraint limits,
% the constraint now becoming:
%
%   xN_min - x <= v_i
%   x - xN_max <= v_i
% 
% Even if small violations on the soft constraints can be tolerated,
% sometimes there are phisical limits which cannot be surpased. In those 
% cases we can enable hard limits such that:
%
% xN_min_hard <= xN_min <= x <= xN_max <= xN_max_hard
%
% In:
%   - mpc: CHRONOS mpc structure
%
%   - x_ter_min (optional): nx column vector, lower bound constraint values
%   on the terminal state vector
%
%   - x_ter_max (optional): nx column vector, upper bound constraint values
%   on the terminal state vector
%
%   - x_ter_min_slack_active (optional): single boolean or nx boolean column 
%   vector, indicates which elements of the state vector minumum constraints
%   have slack variables
%
%   - x_ter_max_slack_active (optional): single boolean or nx boolean column 
%   vector, indicates which elements of the state vector maximum constraints
%   have slack variables
%
%   - x_ter_min_hard (optional): nx column vector, maximum lower bound 
%   constraint values on the state vector
%
%   - x_ter_max_hard (optional): nx column vector, maximum upper bound 
%   constraint values on the state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max,...
                           x_ter_min_slack_active,x_ter_max_slack_active,...
                           x_ter_min_hard,x_ter_max_hard)
arguments
    mpc
    x_ter_min
    x_ter_max
    x_ter_min_slack_active = []
    x_ter_max_slack_active = []
    x_ter_min_hard = []
    x_ter_max_hard = []
end

s_ter_cnstr.min = x_ter_min; 
s_ter_cnstr.max = x_ter_max;

is_there_slack = 0;

if ~isempty(s_ter_cnstr.min)

    s_ter_cnstr.fi_min_x0 = zeros(mpc.nx,1);

    s_ter_cnstr.grad_min = -1 * genGradXter(mpc.N,mpc.N_ctr_hor,...
    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);
    
    % consider slack variable on the gradient
    if ~isempty(x_ter_min_slack_active)

        is_there_slack = 1;

        [mpc,s_ter_cnstr] = init_slack_min_condition(mpc,s_ter_cnstr,...
                            x_ter_min_slack_active,...
                            x_ter_min_hard,mpc.nx,mpc.nx);
    else
        s_ter_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [s_ter_cnstr.hess_min,mi] = genHessIneq(s_ter_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    s_ter_cnstr.grad_min_feas_slv = [s_ter_cnstr.grad_min;-ones(1,mpc.nx)];
    
    s_ter_cnstr.hess_min_feas_slv = genHessIneq(s_ter_cnstr.grad_min_feas_slv);

end

if ~isempty(s_ter_cnstr.max)

    s_ter_cnstr.fi_max_x0 = zeros(mpc.nx,1);

    s_ter_cnstr.grad_max = genGradXter(mpc.N,mpc.N_ctr_hor,...
                            mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);
    
    % consider slack variable on the gradient
    if ~isempty(x_ter_max_slack_active)

        is_there_slack = 1;

        [mpc,s_ter_cnstr] = init_slack_max_condition(mpc,s_ter_cnstr,...
                            x_ter_max_slack_active,...
                            x_ter_max_hard,mpc.nx,mpc.nx);
    else
        s_ter_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [s_ter_cnstr.hess_max,mi] = genHessIneq(s_ter_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    s_ter_cnstr.grad_max_feas_slv = [s_ter_cnstr.grad_max;-ones(1,mpc.nx)];
    
    s_ter_cnstr.hess_max_feas_slv = genHessIneq(s_ter_cnstr.grad_max_feas_slv);
    
end

mpc.s_ter_cnstr = s_ter_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end