function [mpc,feas] = get_terConst_slack_update(mpc,x_ref)

% compute error at terminal state
x_err = x_ref-mpc.s_ter;

fi_ter = x_err'*mpc.P*x_err-1;

mpc.v_ter = max([mpc.slack_ter_epsilon*2,...
                 fi_ter + 1]);

mpc.v(mpc.v_ter_global_index) = mpc.v_ter;

% Terminal Set Constraint Evaluation
mpc.fi_ter_x0 = fi_ter-mpc.v_ter;

mpc.fi_ter_slack_positivity_x0 = mpc.slack_ter_epsilon - mpc.v_ter;

feas = mpc.fi_ter_x0<0 && mpc.fi_ter_slack_positivity_x0<0; 
end