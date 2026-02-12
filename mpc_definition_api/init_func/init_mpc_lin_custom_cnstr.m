%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max,...
%           h_min_slack_active,h_max_slack_active,h_min_hard,h_max_hard)
%
% Define constraints on a custom user defined signal h. The signal h is
% defined as:
%
%   h = Ch * x + Dh * u + Ddh * dh
%
% The user must define the signal h by selecting appropiate values for the
% matrices Ch, Dh, Ddh.
% 
% The constraints on h are then of the form:
%
%   h_min <= h <= h_max
%
% In some scenarios, it is ok if these limits are violated by a small
% amount. In these cases we can enable slack variables v_i on the
% constraints to allow for violations of the now soft constraint limits,
% the constraint now becoming:
%
%   h_min - h <= v_i
%   h - h_max <= v_i
% 
% Even if small violations on the soft constraints can be tolerated,
% sometimes there are phisical limits which cannot be surpased. In those 
% cases we can enable hard limits such that:
%
% h_min_hard <= h_min <= x <= h_max <= h_max_hard
%
% In:
%   - mpc: CHRONOS mpc structure
%
%   - Ch: nh x nx matrix, states output matrix
%
%   - Dh: nh x nu matrix, input feedtrhough matrix
%
%   - Ddh: nh x ndh matrix, disturbance feedtrhough matrix
%
%   - h_min: nh column vector, lower bound constraint values on the user
%   defined signal h
%
%   - h_max: nh column vector, upper bound constraint values on the user
%   defined signal h
%
%   - h_min_slack_active (optional): single boolean or nh boolean column 
%   vector, indicates which elements of the user defined signal minumum 
%   constraints have slack variables
%
%   - h_max_slack_active (optional): single boolean or nh boolean column 
%   vector, indicates which elements of the user defined signal maximum 
%   constraints have slack variables
%
%   - h_min_hard (optional): nh column vector, maximum lower bound 
%   constraint values on the user defined signal
%
%   - h_max_hard (optional): nh column vector, maximum upper bound 
%   constraint values on the user defined signal
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
% Example:
% We want to limit the variation of the MPC control action with respect a
% a given value u_star, e.g. we want to constraint h = u - u_star.
% To achieve this we define h by selecting:
%   - Ch = [0 0 ... 0]
%   - Dh = [1]
%   - Ddh = [-1]
% which corresponds to h = [0 0 ... 0] * x + [1] * u + [-1] * dh
% The "disturbance" term on h corresponds to u_star, to be introduced on 
% the appropiate field on mpc_solve() during runtime MPC execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max, ...
                h_min_slack_active,h_max_slack_active,h_min_hard,h_max_hard)
arguments
    mpc
    Ch
    Dh
    Ddh
    h_min
    h_max
    h_min_slack_active = []
    h_max_slack_active = []
    h_min_hard = []
    h_max_hard = []
end

% General Inequality Matrix
mpc.Ch = Ch;
mpc.Dh = Dh;
mpc.Ddh = Ddh;

mpc.ndh = size(Ddh,2);  %number of disturbance inputs to general inequalities
mpc.nh = size(Ch,1);  %number of general inequalities

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
                        h_min_slack_active,h_min_hard,mpc.Nh,mpc.nh);
    else
        h_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [h_cnstr.hess_min,mi] = genHessIneq(h_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    h_cnstr.grad_min_feas_slv = [h_cnstr.grad_min;-ones(1,mpc.Nh)];
    
    h_cnstr.hess_min_feas_slv = genHessIneq(h_cnstr.grad_min_feas_slv);

end
if ~isempty(h_cnstr.max)

    h_cnstr.fi_max_x0 = zeros(mpc.Nh,1);

    h_cnstr.grad_max = genGradY(mpc.Ch,mpc.Dh,mpc.N,mpc.N_ctr_hor,...
                       mpc.Nx,mpc.Nu,mpc.Nh,mpc.nx,mpc.nu,mpc.nh,mpc.Nv);

    % consider slack variable on the gradient
    if ~isempty(h_max_slack_active)

        is_there_slack = 1;

        [mpc,h_cnstr] = init_slack_max_condition(mpc,h_cnstr,...
                        h_max_slack_active,h_max_hard,mpc.Nh,mpc.nh);
    else
        h_cnstr.min_slack_nv = 0;
    end

    % hessian created after slack is considered on the gradient
    [h_cnstr.hess_max,mi] = genHessIneq(h_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    h_cnstr.grad_max_feas_slv = [h_cnstr.grad_max;-ones(1,mpc.Nh)];
    
    h_cnstr.hess_max_feas_slv = genHessIneq(h_cnstr.grad_max_feas_slv);

end

mpc.h_cnstr = h_cnstr;

if is_there_slack
    mpc = expand_gradients_hessians(mpc);
end

end