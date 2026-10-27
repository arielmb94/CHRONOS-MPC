function mpc = get_mpc_x(mpc,packed_primal)

    mpc.x(:,:)= packed_primal(mpc.x_index_k);
    mpc.x_ter(:) = packed_primal(mpc.x_index_k(:,mpc.N));
        
end