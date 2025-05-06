function x = arrange_opt_vec(s,u,mpc)

nu = mpc.nu;
nx = mpc.nx;

if length(s)<mpc.Nx
    s = fill_vec(s,nu,mpc.Nx,0);
end
if length(u)<mpc.Nx
    u = fill_vec(u,nu,mpc.Nu,0);
end

x = zeros(mpc.Nu+mpc.Nx,1);
x(1:nu) = u(1:nu);
for k = 1:mpc.N
    if k < mpc.N_ctr_hor
        s_k = s(1+(k-1)*nx:k*nx);
        u_k = u(nu+1+(k-1)*nu:nu+k*nu);
        x(nu+1+(nu+nx)*(k-1):nu+(nu+nx)*k) = [s_k;u_k];
    else
        s_k = s(1+(k-1)*nx:k*nx);
        x(mpc.Nu+(mpc.N_ctr_hor-1)*nx+1+nx*(k-mpc.N_ctr_hor):...
            mpc.Nu+(mpc.N_ctr_hor-1)*nx+nx+nx*(k-mpc.N_ctr_hor)) = s_k;
    end
end

end

