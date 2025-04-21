function [s,s_all,s_ter] = get_x(x,s_prev,nx,nu,N,N_h_ctr,Nx)

    s = zeros(Nx,1);
    for k = 1:N
        if k < N_h_ctr
            s((k-1)*nx+1:k*nx) = x(nx*(k-1) + nu*k + 1 : nx*k + nu*k);
        else
            s((k-1)*nx+1:k*nx) = x(nu+(nu+nx)*(N_h_ctr-1)+nx*(k-N_h_ctr)+1 : ...
                nu+(nu+nx)*(N_h_ctr-1)+nx*(k-N_h_ctr)+nx);
        end
    
    end

    s_all = [s_prev;s];
    s_ter = x(end-nx+1 : end);
        
end