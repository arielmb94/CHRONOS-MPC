function du = get_mpc_diff_u(du,u,su,u_prev,N)

    du(:,1) = u(:,1)-u_prev;
    % delta u needs to be computed with su to match gradient definition
    du(:,2:N) = u(:,2:N)-su(:,1:N-1);

end
