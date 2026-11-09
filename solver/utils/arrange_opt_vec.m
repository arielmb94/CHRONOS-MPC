function packed_primal = arrange_opt_vec(x,u,mpc)

nu = mpc.nu;
nx = mpc.nx;

if length(x)<mpc.Nx
    x = fill_vec(x,nu,mpc.Nx,0);
end
if length(u)<mpc.Nx
    u = fill_vec(u,nu,mpc.Nu,0);
end

packed_primal = zeros(mpc.Nu+mpc.Nx,1);
packed_primal(1:nu) = u(1:nu);
for k = 1:mpc.N
    if k < mpc.N_ctr_hor
        x_k = x(1+(k-1)*nx:k*nx);
        u_k = u(nu+1+(k-1)*nu:nu+k*nu);
        packed_primal(nu+1+(nu+nx)*(k-1):nu+(nu+nx)*k) = [x_k;u_k];
    else
        x_k = x(1+(k-1)*nx:k*nx);
        packed_primal(mpc.Nu+(mpc.N_ctr_hor-1)*nx+1+nx*(k-mpc.N_ctr_hor):...
            mpc.Nu+(mpc.N_ctr_hor-1)*nx+nx+nx*(k-mpc.N_ctr_hor)) = x_k;
    end
end

end

