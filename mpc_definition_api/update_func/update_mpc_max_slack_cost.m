% UPDATE_MPC_MAX_SLACK_COST Updates the penalty weights for slack values on the 
% maximum limits of box constraints with soft constraints defined during
% initialization.
%
%   mpc = UPDATE_MPC_MAX_SLACK_COST(mpc, cnstr, qv_max) sets to qv_max the
%   penalty weights of the slack values active on the maximum limits of the
%   given constraint and adapts the CHRONOS mpc gradients accordignly. 
%
%   The list of available CHRONOS box constraint structures is:
%
%   mpc.s_cnstr            - box constraint on the states
%   mpc.s_ter_cnstr        - box constraint on the terminal state
%   mpc.u_cnstr            - box constraint on the control action
%   mpc.du_cnstr           - box constraint on the control action rate
%   mpc.y_cnstr            - box constraint on the output tracking signals
%   mpc.h_cnstr            - box constraint on the user defined constraints
% 
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       cnstr              - CHRONOS structure for the box constraint to be
%                          updated.                         
%       qv_min             - [ni x 1] or scalar. Penalty weight for 
%                          violating the cnstr maximum soft limits. Higher 
%                          values mean stricter enforcement. ni is the full
%                          size of the constraint (e.g. nx for constraints
%                          on the state vector, nu for constraints on the
%                          control actions and control action rate, ny for
%                          constraints of the tracking signal and nh for
%                          user defined contraints)
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary 
%                          background math are automatically assembled and 
%                          added to the object.
%
%  USAGE TIPS:
%       - To be used only on constraints that already have soft constraints
%       defined during initialization of the CRHONOS mpc problem
%       - Passing a scalar to qv_max will automatically apply that setting 
%       across all slacks active on the soft constraint.
%       - Passing an [ni x 1] vector you can modify individually the penalty
%       weight for each constraint element.
function mpc = update_mpc_max_slack_cost(mpc,cnstr,qv_max_in)

n_i = length(cnstr.max);
qv_max = zeros(n_i,1);
if length(qv_max_in)==1 && n_i > 1
    qv_max(:) = qv_max_in*ones(n_i,1);
else
    qv_max(:) = qv_max_in;
end

mpc.gradSlackqv(mpc.Nx+mpc.Nu+cnstr.max_v_global_index) = ...
        cnstr.max_slack_local_map*qv_max;

end