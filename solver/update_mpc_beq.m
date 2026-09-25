function mpc = update_mpc_beq(mpc,s_prev,A,dyn_use_d)

mpc.beq_0(:) = -A*s_prev;

if ~isempty(dyn_use_d)
    mpc = beq_disturbance_wrapper(mpc);
end
end
