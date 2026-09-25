function mpc = beq_disturbance_wrapper(mpc)
[mpc.beq_0,mpc.beq_k] = beq_disturbance_local(mpc.beq_0,mpc.beq_k,mpc.Bd,mpc.d,mpc.N);
end

function [beq_0,beq_k] = beq_disturbance_local(beq_0,beq_k,Bd,d,N)
beq_0(:) = beq_0(:) - Bd(:,:,1)*d(:,1);
for k = 1:N-1
    beq_k(:,k) = -Bd(:,:,k+1)*d(:,k+1);
end
end
