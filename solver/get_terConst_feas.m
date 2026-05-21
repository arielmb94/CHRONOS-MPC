function [mpc,feas] = get_terConst_feas(mpc,x_ref)

% compute error at terminal state
x_err = x_ref-mpc.s_ter;

% Terminal Set Constraint Evaluation
mpc.fi_ter_x0 = x_err'*mpc.P*x_err-1-mpc.v_ter;

mpc.fi_ter_slack_positivity_x0 = mpc.slack_ter_epsilon - mpc.v_ter;
tol = -1e-4;
feas = mpc.fi_ter_x0<tol && mpc.fi_ter_slack_positivity_x0<0; 
end