function mpc = update_mpc_f0_hess(mpc,tracking_cost,...
    quad_control_cost,controlrate_cost,quad_custom_cost,ter_ingredients)

    mpc = f0_hess_reset_wrapper(mpc);
    if ~isempty(tracking_cost), mpc = f0_hess_tracking_wrapper(mpc); end
    if ~isempty(quad_control_cost), mpc = f0_hess_control_wrapper(mpc); end
    if ~isempty(controlrate_cost), mpc = f0_hess_control_rate_wrapper(mpc); end
    if ~isempty(quad_custom_cost), mpc = f0_hess_custom_quad_wrapper(mpc); end
    if ~isempty(ter_ingredients), mpc = f0_hess_terminal_wrapper(mpc); end

end
