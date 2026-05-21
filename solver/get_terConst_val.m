function mpc = get_terConst_val(mpc,x_ref,s_ter,P,v_feas)

% compute error at terminal state
x_err = x_ref-s_ter;

% Terminal Set Constraint Evaluation
fi_ter_x0 = x_err'*P*x_err;

% compute slack variable value
mpc.v_ter = max(1 + mpc.slack_ter_epsilon, fi_ter_x0 + mpc.slack_ter_epsilon);

% Map to global v vector
mpc.v = mpc.v + mpc.ter_cnstr_map * mpc.v_ter;

% 2. Apply the updated slack variables across all time instances
mpc.fi_ter_x0 = fi_ter_x0 - mpc.v_ter;

% for feasibility solver: fi - v_feas <= 0
if v_feas
    fi_ter_x0 = fi_ter_x0 - v_feas;
end

end