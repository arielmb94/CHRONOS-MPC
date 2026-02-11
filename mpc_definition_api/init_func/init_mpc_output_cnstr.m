%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_output_cnstr(mpc,y_min,y_max,...
%         y_min_slack_active,y_max_slack_active,y_min_hard,y_max_hard))
%
% Define constraints on the tracking signal y:
%
%   y_min <= y <= y_max
%
% In some scenarios, it is ok if these limits are violated by a small
% amount. In these cases we can enable slack variables v_i on the
% constraints to allow for violations of the now soft constraint limits,
% the constraint now becoming:
%
%   y_min - y <= v_i
%   y - y_max <= v_i
% 
% Even if small violations on the soft constraints can be tolerated,
% sometimes there are phisical limits which cannot be surpased. In those 
% cases we can enable hard limits such that:
%
% y_min_hard <= y_min <= y <= y_max <= y_max_hard
%
% In:
%   - mpc: CHRONOS mpc structure
%
%   - y_min (optional): ny column vector, lower bound constraint values on 
%   the tracking signal
%
%   - y_max (optional): ny column vector, upper bound constraint values on 
%   the tracking signal
%
%   - y_min_slack_active (optional): single boolean or ny boolean column 
%   vector, indicates which elements of the output vector minumum constraints
%   have slack variables
%
%   - y_max_slack_active (optional): single boolean or ny boolean column 
%   vector, indicates which elements of the output vector maximum constraints
%   have slack variables
%
%   - y_min_hard (optional): ny column vector, maximum lower bound 
%   constraint values on the state vector
%
%   - y_max_hard (optional): nx column vector, maximum upper bound 
%   constraint values on the state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_output_cnstr(mpc,y_min,y_max,...
                y_min_slack_active,y_max_slack_active,y_min_hard,y_max_hard)
arguments
    mpc
    y_min
    y_max
    y_min_slack_active = []
    y_max_slack_active = []
    y_min_hard = []
    y_max_hard = []
end

y_cnstr.min = y_min;
y_cnstr.max = y_max;

is_there_slack = 0;

if ~isempty(y_cnstr.min)

    y_cnstr.fi_min_x0 = zeros(mpc.Ny,1);

    % Outputs box constraints
    y_cnstr.grad_min = -1 * genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
                            mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

    % consider slack variable on the gradient
    if ~isempty(y_min_slack_active)

        is_there_slack = 1;

        [mpc,y_cnstr] = init_slack_min_condition(mpc,y_cnstr,...
                        y_min_slack_active,y_min_hard,mpc.Ny,mpc.ny);
    else
        y_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [y_cnstr.hess_min,mi] = genHessIneq(y_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    y_cnstr.grad_min_feas_slv = [y_cnstr.grad_min;-ones(1,mpc.Ny)];
    
    y_cnstr.hess_min_feas_slv = genHessIneq(y_cnstr.grad_min_feas_slv);

end

if ~isempty(y_cnstr.max)

    y_cnstr.fi_max_x0 = zeros(mpc.Ny,1);

    % Outputs box constraints
    y_cnstr.grad_max = genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
                       mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

    % consider slack variable on the gradient
    if ~isempty(y_max_slack_active)

        is_there_slack = 1;

        [mpc,y_cnstr] = init_slack_max_condition(mpc,y_cnstr,...
                        y_max_slack_active,y_max_hard,mpc.Ny,mpc.ny);
    else
        y_cnstr.max_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [y_cnstr.hess_max,mi] = genHessIneq(y_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    y_cnstr.grad_max_feas_slv = [y_cnstr.grad_max;-ones(1,mpc.Ny)];
    
    y_cnstr.hess_max_feas_slv = genHessIneq(y_cnstr.grad_max_feas_slv);

end

mpc.y_cnstr = y_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end