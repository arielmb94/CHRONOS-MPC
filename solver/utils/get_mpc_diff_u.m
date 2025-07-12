function mpc = get_mpc_diff_u(u_prev,mpc)

    u_total = [u_prev reshape(mpc.u,[mpc.nu,mpc.N_ctr_hor])];
    
    delta_u = diff(u_total,1,2);

    mpc.du(:) = reshape(delta_u,[mpc.Nu 1]);
        
end