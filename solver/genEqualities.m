function mpc = genEqualities(mpc,A,B,N,N_h_ctr,nx,nu)

for k = 0:N-1

    if k == 0
        mpc.Aeq(1:nx,1:nu+nx) = [B -eye(nx)];
    
    elseif k < N_h_ctr
        mpc.Aeq(k*nx+1:(k+1)*nx,(nu+nx)*k+1-nx:(nu+nx)*(k+1)) = [A B -eye(nx)];

    else
        mpc.Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+1:(nx+nu)*(N_h_ctr-1)+nu) = ...
            B;

        mpc.Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr))+1:...
            (nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr-2))) = ...
            [A -eye(nx)];

    end
end

end