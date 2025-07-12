function mpc = get_mpc_x(x,s_prev,mpc)

    for k = 1:mpc.N
        if k < mpc.N_ctr_hor
            mpc.s((k-1)*mpc.nx+1:k*mpc.nx) = x(mpc.nx*(k-1) + mpc.nu*k + 1 : mpc.nx*k + mpc.nu*k);
        else
            mpc.s((k-1)*mpc.nx+1:k*mpc.nx) = x(mpc.nu+(mpc.nu+mpc.nx)*(mpc.N_ctr_hor-1)+mpc.nx*(k-mpc.N_ctr_hor)+1 : ...
                mpc.nu+(mpc.nu+mpc.nx)*(mpc.N_ctr_hor-1)+mpc.nx*(k-mpc.N_ctr_hor)+mpc.nx);
        end
    
    end

    mpc.s_all = [s_prev;mpc.s];
    mpc.s_ter(:) = mpc.s((mpc.Nx)-mpc.nx+1 : mpc.Nx);
        
end