% INIT_MPC_TER_INGREDIENTS_DLQR Initializes terminal cost and constraints for MPC
%
% This function computes the terminal ingredients required to guarantee closed-loop 
% stability in the MPC formulation. It uses the Discrete Linear Quadratic Regulator 
% (DLQR) to calculate the infinite-horizon cost-to-go matrix (P) and the optimal 
% terminal feedback gain (K). The terminal cost adds the following penalty on the
% final state of the MPC Prediction Horizon:
%
%   J += (ref_xN - xN)' * P * (ref_xN - xN)
% 
% (Optional, but not recommended) It can also configure a terminal set constraint:
%
%   (ref_xN - xN)' * P * (ref_xN - xN) < 1
%
% Usage:
%   mpc = init_mpc_ter_ingredients_dlqr(mpc, Qx, Ru, x_ref_is_y, terminal_constraint, qv_ter)
%
% IMPORTANT: before calling this function you must ensure that the dynamics
% of the mpc have been initialized with representative values. Proper
% initialization can be done during the call to init_mpc_system() or by 
% modifying the MPC struct fields mpc.A and mpc.B manually.
%
% Inputs:
%   mpc                 - CHRONOS mpc structure.
%   Qx                  - [nx x nx] State weighting matrix for the DLQR calculation.
%   Ru                  - [nu x nu] Control weighting matrix for the DLQR calculation.
%   x_ref_is_y          - (Optional) Boolean flag. In cases where the tracking signal 
%                         is set to thefull state vector, e.g. y = I * x, setting 
%                         x_ref_is_y to 1 allows you to avoid defining ref_xN as a 
%                         seperate argument on mpc_solve() during runtime MPC calls.
%                         Default is 0.
%   terminal_constraint - (Optional) Boolean flag to enable the terminal set constraint.
%                         Default is 0 (Disabled).
%   qv_ter              - (Optional) Penalty weight for the terminal constraint slack 
%                         variable. If left empty, it defaults to the global mpc.qv.
%
% Outputs:
%   mpc                 - Updated MPC structure. All necessary background math 
%                         (constraint gradients, Hessians, and slack variables) 
%                         are automatically assembled and added to the object.
%
% ==============================================================================
% ⚠️PRACTICAL USAGE NOTE: Terminal Constraints
% ==============================================================================
% The 'terminal_constraint' option is included primarily for theoretical 
% completeness (e.g., rigid proofs of recursive feasibility and stability). 
% 
% In practice, ENABLING THIS CONSTRAINT IS NOT RECOMMENDED. If the terminal 
% region of attraction is too small, it forces the solver against a steep 
% log-barrier, causing severe numerical stiffness, Hessian blowups, and erratic 
% control behavior. 
%
% For the vast majority of real-world applications, the terminal cost alone 
% provides more than enough "pull" to stabilize the system and guide it 
% without sacrificing solver speed or numerical stability.
function mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,...
                                             x_ref_is_y,...
                                             terminal_constraint,qv_ter)
arguments
    mpc
    Qx
    Ru
    x_ref_is_y = 0
    terminal_constraint = 0
    qv_ter = []
end

mpc.ter_ingredients = 1;
mpc.ter_constraint = terminal_constraint;
mpc.x_ref_is_y = x_ref_is_y;

[K,P] = dlqr(mpc.A,mpc.B,Qx,Ru);

mpc.K = K;
mpc.P = P;
mpc.P2 = 2*P;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx+mpc.Nv);
end

mpc.hessTerminalCost = zeros(mpc.Nx+mpc.Nu+mpc.Nv);
N = mpc.Nx+mpc.Nu;
mpc.hessTerminalCost(N-mpc.nx+1: N,N-mpc.nx+1 : N) = 2*P;

mpc.hessCost = mpc.hessCost + mpc.hessTerminalCost;

if terminal_constraint 
    % fi_ter = e_xN'*P*e_xN - vi < 0
    % grad_fi = 2*grad_e_xN*P*e_xN - [0;1]
    % hess_fi = 2*grad_e_xN*P*grad_e_xN' (same as terminal cost function)

    mpc.fi_ter_x0 = 0;

    mpc.v = [mpc.v;0];
    mpc.v_ter = 0;
    % get index of v_ter in global v vector
    mpc.v_ter_global_index = mpc.Nv+1;

    % gradient/hess of constraint: -v<=0 (slack positivity constraint)
    mpc.fi_ter_slack_positivity_x0 = 0;
    mpc.fi_ter_slack_positivity_grad = -1 * [zeros(mpc.Nx+mpc.Nu+mpc.Nv,1);...
                                            1];
    [mpc.fi_ter_slack_positivity_hess,mi] = genHessIneq(mpc.fi_ter_slack_positivity_grad);
    mpc.m = mpc.m+mi;

    % Initialize Penalty term for new slack variables
    if isempty(qv_ter)
        qv_ter = mpc.qv;
    end

    if ~isempty(mpc.gradSlackqv)
        % expand slack penalty term grad
        mpc.gradSlackqv = [mpc.gradSlackqv;
                            qv_ter];
    else
        mpc.gradSlackqv = [zeros(mpc.Nx+mpc.Nu,1);
                            qv_ter];
    end

    % update global counter of slack variables
    mpc.Nv = mpc.Nv+1;
    % adapt gradients/hessians for the new variable vector size
    mpc = expand_gradients_hessians(mpc);
end

end