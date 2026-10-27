function du = get_mpc_diff_u(du,u,xu,u_prev,N)

    du(:,1) = u(:,1)-u_prev;
    % delta u needs to be computed with xu to match gradient definition
    du(:,2:N) = u(:,2:N)-xu(:,1:N-1);

end
