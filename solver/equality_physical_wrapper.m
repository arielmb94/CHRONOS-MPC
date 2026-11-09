function mpc = equality_physical_wrapper(mpc)
[mpc.rp_0,mpc.rp_k] = equality_physical_local(mpc.rp_0,mpc.rp_k,mpc.A,mpc.B,...
    mpc.beq_0,mpc.beq_k,mpc.u,mpc.x,mpc.x_col,mpc.N);
end

function [rp_0,rp_k] = equality_physical_local(rp_0,rp_k,A,B,beq_0,beq_k,u,x,x_idx,N)
% k = 0
rp_0(x_idx) = B(:,:,1)*u(:,1)-x(:,1)-beq_0;

for k = 1:N-1
    % [A B -I][x u x+]'-beq
    rp_k(x_idx,k) = A(:,:,k+1)*x(:,k)+B(:,:,k+1)*u(:,k+1)-x(:,k+1)...
                            -beq_k(:,k);
end
end
