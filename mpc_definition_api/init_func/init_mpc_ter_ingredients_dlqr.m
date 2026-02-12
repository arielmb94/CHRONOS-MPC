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
                                             terminal_constraint,x_ref_is_y)

mpc.ter_ingredients = 1;
mpc.ter_constraint = terminal_constraint;
mpc.x_ref_is_y = x_ref_is_y;

[~,P] = dlqr(mpc.A,mpc.B,Qx,Ru);

mpc.P = P;
mpc.P2 = 2*P;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx+mpc.Nv);
end

mpc.fi_ter_x0 = 0;
mpc.hessTerminalCost = zeros(mpc.Nx+mpc.Nu+mpc.Nv);
N = mpc.Nx+mpc.Nu;
mpc.hessTerminalCost(N-mpc.nx+1: N,N-mpc.nx+1 : N) = 2*P;

mpc.hessCost = mpc.hessCost + mpc.hessTerminalCost;

if terminal_constraint 
    % fi_ter = e_xN'*P*e_xN - vi < 0
    % grad_fi = 2*grad_e_xN*P*e_xN - [0;1]
    % hess_fi = 2*grad_e_xN*P*grad_e_xN' (same as terminal cost function)

    mpc.ter_cnstr_slack_index = [zeros(1,mpc.Nx+mpc.Nu+mpc.Nv) 1];

    % initialize terminal constraint gradient with constant slack element
    % -> -[0;1]
    mpc.ter_cnstr_grad = [zeros(mpc.Nx+mpc.Nu+mpc.Nv,1);
                            -1];

    % Initialize Penalty term for new slack variables -> Jv = vi'*Qv*vi
    if ~isempty(mpc.gradSlackQv)
        % expand slack penalty term grad
        mpc.gradSlackQv = [mpc.gradSlackQv;
                           zeros(1,mpc.Nv)];
        % add columns for new slack
        mpc.gradSlackQv = [mpc.gradSlackQv mpc.ter_cnstr_slack_index'*mpc.Qv];

        % expand slack penalty term hess
        mpc.hessSlackTerm = [mpc.hessSlackTerm zeros(mpc.Nx+mpc.Nu+mpc.Nv,1);
            zeros(1,mpc.Nx+mpc.Nu+mpc.Nv) mpc.Qv];
    else
        mpc.gradSlackQv = cnstr.min_slack_index'*mpc.Qv;
        mpc.hessSlackTerm = [zeros(mpc.Nx+mpc.Nu) zeros(mpc.Nx+mpc.Nu,1);
            zeros(1,mpc.Nx+mpc.Nu) mpc.Qv];
    end

    % update global counter of slack variables
    mpc.Nv = mpc.Nv+1;
    % adapt gradients/hessians for the new variable vector size
    mpc = expand_gradients_hessians(mpc);
end

end