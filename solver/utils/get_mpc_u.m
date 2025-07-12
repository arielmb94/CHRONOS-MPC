function mpc = get_mpc_u(x,mpc)

    for k = 1:mpc.N_ctr_hor
        
        mpc.u((k-1)*mpc.nu+1:k*mpc.nu) = x(mpc.nu + mpc.nx*(k-1) + mpc.nu*(k-2)+1:...
                                            mpc.nu + mpc.nx*(k-1) + mpc.nu*(k-1));
    
    end
        
end