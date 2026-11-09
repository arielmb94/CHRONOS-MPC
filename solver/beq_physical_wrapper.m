function mpc = beq_physical_wrapper(mpc,x_prev)
mpc.beq_0 = beq_physical_local(mpc.beq_0,mpc.A,x_prev);
end

function beq_0 = beq_physical_local(beq_0,A,x_prev)
%% dynamics
% k = 0
beq_0(:) = -A(:,:,1)*x_prev;
end
