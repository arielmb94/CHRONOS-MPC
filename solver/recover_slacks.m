function mpc = recover_slacks(mpc,has_s_cnstr,has_u_cnstr,has_du_cnstr,...
    has_y_cnstr,has_h_cnstr)
if ~isempty(has_s_cnstr)
    mpc = recover_slacks_state_wrapper(mpc,mpc.s_cnstr);
end
if ~isempty(has_u_cnstr)
    mpc = recover_slacks_control_wrapper(mpc,mpc.u_cnstr);
end
if ~isempty(has_du_cnstr)
    mpc = recover_slacks_control_rate_wrapper(mpc,mpc.du_cnstr);
end
if ~isempty(has_y_cnstr)
    mpc = recover_slacks_output_wrapper(mpc,mpc.y_cnstr);
end
if ~isempty(has_h_cnstr)
    mpc = recover_slacks_custom_wrapper(mpc,mpc.h_cnstr);
end
end
