function cnstr = init_min_empty(cnstr)

cnstr.min_limit = 0;
cnstr.fi_min_x0 = [];
cnstr.grad_min = [];
cnstr.hess_min = genHessIneq(cnstr.grad_min);
cnstr.min_v = [];
cnstr.min_v_global_index = [];
cnstr.min_slack_positivity_fi_x0 = [];
cnstr.min_slack_positivity_grad = [];
cnstr.min_slack_positivity_hess = genHessIneq(cnstr.min_slack_positivity_grad);
end