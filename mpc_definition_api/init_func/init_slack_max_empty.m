function cnstr = init_slack_max_empty(cnstr)

cnstr.max_slack_nv = 0;
cnstr.max_slack_active = [];
cnstr.max_v = [];
cnstr.max_slack_local_map = [];
cnstr.max_v_global_index = [];
cnstr.max_slack_positivity_fi_x0 = [];
cnstr.max_slack_positivity_grad = [];
cnstr.max_slack_positivity_hess = genHessIneq(cnstr.max_slack_positivity_grad);
end