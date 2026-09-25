function mpc = grad_f0_MPC(mpc,quad_control_cost,lin_control_cost,controlrate_cost,...
    tracking_cost,quad_custom_cost,lin_custom_cost,ter_ingredients)

mpc = grad_f0_reset_wrapper(mpc);
if ~isempty(quad_control_cost) || ~isempty(lin_control_cost), mpc = grad_f0_control_wrapper(mpc,quad_control_cost,lin_control_cost); end
if ~isempty(controlrate_cost), mpc = grad_f0_control_rate_wrapper(mpc); end
if ~isempty(tracking_cost), mpc = grad_f0_tracking_wrapper(mpc); end
if ~isempty(quad_custom_cost), mpc = grad_f0_custom_quad_wrapper(mpc); end
if ~isempty(lin_custom_cost), mpc = grad_f0_custom_lin_wrapper(mpc); end
if ~isempty(ter_ingredients), mpc = grad_f0_terminal_wrapper(mpc); end
mpc = grad_f0_scale_wrapper(mpc);
end
