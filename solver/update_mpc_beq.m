function mpc = update_mpc_beq(mpc,x_prev,A,dyn_use_d)

mpc.beq_0(:) = -A*x_prev;

if ~isempty(dyn_use_d)
    mpc = beq_disturbance_wrapper(mpc);
end
end
