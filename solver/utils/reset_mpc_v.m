function mpc = reset_mpc_v(mpc)
    
    % set to 0 global slack variable vector
    mpc.v(:) = 0;

end