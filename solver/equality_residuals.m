function mpc = equality_residuals(mpc,has_du)

[mpc.rp_0,mpc.rp_k] = equality_physical_local(mpc.rp_0,mpc.rp_k,mpc.A,mpc.B,...
    mpc.beq_0,mpc.beq_k,mpc.u,mpc.s,mpc.s_col,mpc.N);

if ~isempty(has_du)
    [mpc.rp_0,mpc.rp_k] = equality_du_local(mpc.rp_0,mpc.rp_k,mpc.u,mpc.su,...
    mpc.su_col,mpc.N);
end

end

function [rp_0,rp_k] = equality_physical_local(rp_0,rp_k,A,B,beq_0,beq_k,u,s,s_idx,N)
% k = 0
rp_0(s_idx) = B(:,:,1)*u(:,1)-s(:,1)-beq_0;

for k = 1:N-1
    % [A B -I][s u s+]'-beq
    rp_k(s_idx,k) = A(:,:,k+1)*s(:,k)+B(:,:,k+1)*u(:,k+1)-s(:,k+1)...
                            -beq_k(:,k);
end
end

function [rp_0,rp_k] = equality_du_local(rp_0,rp_k,u,su,su_idx,N)
% this is the control input delay equality condition: su+ = u_prev
% [I -I][u su+]' = 0
% k = 0
rp_0(su_idx) = u(:,1)-su(:,1);
% k = 1...N
rp_k(su_idx,:) = u(:,2:N)-su(:,2:N);
end