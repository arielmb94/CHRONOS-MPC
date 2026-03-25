function cnstr = init_slack_min_empty(cnstr)

cnstr.min_slack_nv = 0;
cnstr.min_slack_active = [];
cnstr.min_v = [];
cnstr.min_slack_local_map = [];
cnstr.min_v_global_index = [];
cnstr.min_slack_positivity_fi_x0 = [];
cnstr.min_slack_positivity_grad = [];
cnstr.min_slack_positivity_hess = genHessIneq(cnstr.min_slack_positivity_grad);
end