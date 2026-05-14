function cnstr = init_max_empty(cnstr)

cnstr.max_limit = 0;
cnstr.fi_max_x0 = [];
cnstr.grad_max = [];
cnstr.hess_max = genHessIneq(cnstr.grad_max);
cnstr.max_v = [];
cnstr.max_v_global_index = [];
cnstr.max_slack_positivity_fi_x0 = [];
cnstr.max_slack_positivity_grad = [];
cnstr.max_slack_positivity_hess = genHessIneq(cnstr.max_slack_positivity_grad);
end