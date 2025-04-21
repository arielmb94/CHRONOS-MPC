function x0 = update_mpc_warmstart(x,mpc,fill_style)
arguments
    x
    mpc
    fill_style = 0;
end

    s_prev = zeros(mpc.nx,1); % s_prev is not used in warmstarting
    s = get_x(x,s_prev,mpc.nx,mpc.nu,mpc.N,mpc.N_ctr_hor,mpc.Nx);
    u = get_u(x,mpc.nx,mpc.nu,mpc.N_ctr_hor,mpc.Nu);

    % delete first state step and fill
    s_crop = s(mpc.nx+1:end);
    s_fill = fill_vec(s_crop,mpc.nx,mpc.Nx,fill_style);
    % delete first control step and fill
    u_crop = u(mpc.nu+1:end);
    u_fill = fill_vec(u_crop,mpc.nu,mpc.Nu,fill_style);

    x0 = arrange_opt_vec(s_fill,u_fill,mpc);

end

