function mpc = equality_residuals(mpc,has_du)

[mpc.rp_0,mpc.rp_k] = equality_physical_local(mpc.rp_0,mpc.rp_k,mpc.A,mpc.B,...
    mpc.beq_0,mpc.beq_k,mpc.u,mpc.x,mpc.x_col,mpc.N);

if ~isempty(has_du)
    [mpc.rp_0,mpc.rp_k] = equality_du_local(mpc.rp_0,mpc.rp_k,mpc.u,mpc.xu,...
    mpc.xu_col,mpc.N);
end

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

function [rp_0,rp_k] = equality_du_local(rp_0,rp_k,u,xu,xu_idx,N)
% this is the control input delay equality condition: xu+ = u_prev
% [I -I][u xu+]' = 0
% k = 0
rp_0(xu_idx) = u(:,1)-xu(:,1);
% k = 1...N
rp_k(xu_idx,:) = u(:,2:N)-xu(:,2:N);
end