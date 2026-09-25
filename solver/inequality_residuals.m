function mpc = inequality_residuals(mpc,g_0,g_k,g_ter,has_s_cnstr,has_u_cnstr,...
    has_du_cnstr,has_y_cnstr,has_h_cnstr)

if ~isempty(g_0)
    mpc.ri_0 = inequality_residual_0_local(mpc.ri_0,mpc.g_0,mpc.v_0,mpc.v_rows_0);
end
if ~isempty(g_k)
    mpc.ri_k = inequality_residual_k_local(mpc.ri_k,mpc.g_k,mpc.v_k,mpc.v_rows_k);
end
if ~isempty(g_ter)
    mpc.ri_ter = inequality_residual_terminal_local(mpc.ri_ter,mpc.g_ter,mpc.v_ter);
end
if ~isempty(has_s_cnstr), mpc = inequality_state_wrapper(mpc,mpc.s_cnstr); end
if ~isempty(has_u_cnstr), mpc = inequality_control_wrapper(mpc,mpc.u_cnstr); end
if ~isempty(has_du_cnstr), mpc = inequality_control_rate_wrapper(mpc,mpc.du_cnstr); end
if ~isempty(has_y_cnstr), mpc = inequality_output_wrapper(mpc,mpc.y_cnstr); end
if ~isempty(has_h_cnstr), mpc = inequality_custom_wrapper(mpc,mpc.h_cnstr); end
end

function ri_0 = inequality_residual_0_local(ri_0,g_0,v_0,v_rows_0)
ri_0(:) = g_0;
if ~isempty(v_0)
    ri_0(v_rows_0) = ri_0(v_rows_0) - v_0;
end
end

function ri_k = inequality_residual_k_local(ri_k,g_k,v_k,v_rows_k)
ri_k(:,:) = g_k;
if ~isempty(v_k)
    ri_k(v_rows_k,:) = ri_k(v_rows_k,:) - v_k;
end
end

function ri_ter = inequality_residual_terminal_local(ri_ter,g_ter,v_ter)
ri_ter(:) = g_ter-v_ter;
end
