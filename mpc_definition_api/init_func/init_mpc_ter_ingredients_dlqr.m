%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru, terminal_constraint,x_ref_is_y)
%
% Configures the MPC terminal cost ingredient and optionally enables the
% use of ellipsoidal terminal set contraints.
%
% The terminal cost adds the following penalty on the final state of the 
% MPC Prediction Horizon:
%
%   J += (ref_xN - xN)' * P * (ref_xN - xN)
%
% where the terminal cost weight P is computed internally from the solution
% to a discrete Ricatti equation, obtained from:
%
%   P = dlqr(mpc.A,mpc.B,Qx,Ru),
%
% xN is the final state on the MPC prediction horizon and ref_xN is
% the desired value for the MPC final state. ref_xN is to be introduced on
% mpc_solve() during runtime MPC execution.
%
% If the terminal set constraint is enabled, CHRONOS enforces the following
% ellipsoidal terminal set constraint:
%
%   (ref_xN - xN)' * P * (ref_xN - xN) <= 1
%
% IMPORTANT: before calling this function you must ensure that the dynamics
% of the mpc have been initialized with representative values. Proper
% initialization can be done during the call to init_mpc_system() or by 
% modifying the MPC struct fields mpc.A and mpc.B manually.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Qx: nx x nx matrix, dLQR penalty on the states
%   - Ru: nu x nu matrix, dLQR penalty on the control action
%   - terminal_constraint: Boolean. If set to 1 the terminal set is added
%   to the list of constraints to be satisfied. NOTE: on the current state
%   of CHRONOS it is not advisible to enable terminal set constraints, once
%   slack variables have been added to CHRONOS the feasibility issues will
%   be solved
%   - x_ref_is_y: Boolean. In cases where the tracking signal is set to the
%   full state vector, e.g. y = C * x, setting x_ref_is_y to 1 allows you
%   to avoid defining ref_xN as a seperate argument on mpc_solve() during 
%   runtime MPC calls.
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,...
                                             terminal_constraint,x_ref_is_y, ...
                                             qv_ter)
arguments
    mpc
    Qx
    Ru
    terminal_constraint = 0
    x_ref_is_y = 0
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