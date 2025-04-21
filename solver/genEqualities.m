function Aeq = genEqualities(A,B,N,N_h_ctr,Nx,Nu,nx,nu)

Aeq = zeros(Nx,Nx+Nu);
for k = 0:N-1

    if k == 0
        Aeq(1:nx,1:nu+nx) = [B -eye(nx)];
    
    elseif k < N_h_ctr
        Aeq(k*nx+1:(k+1)*nx,(nu+nx)*k+1-nx:(nu+nx)*(k+1)) = [A B -eye(nx)];

    else
        Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+1:(nx+nu)*(N_h_ctr-1)+nu) = ...
            B;

        Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr))+1:...
            (nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr-2))) = ...
            [A -eye(nx)];

    end
end

end