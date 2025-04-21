function delta_u = diff_u(u,u_prev,nu,N,Nu)

    u_total = [u_prev reshape(u,[nu,N])];
    
    delta_u = diff(u_total,1,2);

    delta_u = reshape(delta_u,[Nu 1]);
        
end