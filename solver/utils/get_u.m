function u = get_u(x,nx,nu,N,Nu)

    u = zeros(Nu,1);
    for k = 1:N
        
        u((k-1)*nu+1:k*nu) = x(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1));
    
    end
        
end