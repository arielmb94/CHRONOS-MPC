function mpc = beq_physical_wrapper(mpc,s_prev)
mpc.beq_0 = beq_physical_local(mpc.beq_0,mpc.A,s_prev);
end

function beq_0 = beq_physical_local(beq_0,A,s_prev)
%% dynamics
% k = 0
beq_0(:) = -A(:,:,1)*s_prev;
end
