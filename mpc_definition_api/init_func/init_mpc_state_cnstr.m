%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_state_cnstr(mpc,x_min,x_max,...
%                              x_min_slack_active,x_max_slack_active,...
%                              x_min_hard,x_max_hard)
%
% Define constraints on the state vector x:
%
%   x_min <= x <= x_max
% 
% In some scenarios, it is ok if these limits are violated by a small
% amount. In these cases we can enable slack variables v_i on the
% constraints to allow for violations of the now soft constraint limits,
% the constraint now becoming:
%
%   x_min - x <= v_i
%   x - x_max <= v_i
% 
% Even if small violations on the soft constraints can be tolerated,
% sometimes there are phisical limits which cannot be surpased. In those 
% cases we can enable hard limits such that:
%
% x_min_hard <= x_min <= x <= x_max <= x_max_hard
%
% In:
%   - mpc: CHRONOS mpc structure
%
%   - x_min (optional): nx column vector, (soft) lower bound constraint values on 
%   the state vector
%
%   - x_max (optional): nx column vector, (soft) upper bound constraint values on 
%   the state vector
%
%   - x_min_slack_active (optional): single boolean or nx boolean column 
%   vector, indicates which elements of the state vector minumum constraints
%   have slack variables
%
%   - x_max_slack_active (optional): single boolean or nx boolean column 
%   vector, indicates which elements of the state vector maximum constraints
%   have slack variables
%
%   - x_min_hard (optional): nx column vector, maximum lower bound 
%   constraint values on the state vector
%
%   - x_max_hard (optional): nx column vector, maximum upper bound 
%   constraint values on the state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_state_cnstr(mpc,x_min,x_max,...
               x_min_slack_active,x_max_slack_active,x_min_hard,x_max_hard)
arguments
    mpc
    x_min
    x_max
    x_min_slack_active = []
    x_max_slack_active = []
    x_min_hard = []
    x_max_hard = []
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
                        x_min_hard,mpc.Nx,mpc.nx);
    else
        s_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [s_cnstr.hess_min,mi] = genHessIneq(s_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    s_cnstr.grad_min_feas_slv = [s_cnstr.grad_min;-ones(1,mpc.Nx)];
    
    s_cnstr.hess_min_feas_slv = genHessIneq(s_cnstr.grad_min_feas_slv);

end

if ~isempty(s_cnstr.max)
    
    s_cnstr.fi_max_x0 = zeros(mpc.Nx,1);

    s_cnstr.grad_max = genGradX(mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(x_max_slack_active)

        is_there_slack = 1;

        [mpc,s_cnstr] = init_slack_max_condition(mpc,s_cnstr,x_max_slack_active,...
                        x_max_hard,mpc.Nx,mpc.nx);
    else
        s_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [s_cnstr.hess_max,mi] = genHessIneq(s_cnstr.grad_max);
    mpc.m = mpc.m+mi;
    
    % initialize feasibility solver max condition
    s_cnstr.grad_max_feas_slv = [s_cnstr.grad_max;-ones(1,mpc.Nx)];
    
    s_cnstr.hess_max_feas_slv = genHessIneq(s_cnstr.grad_max_feas_slv);
    
end

mpc.s_cnstr = s_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end