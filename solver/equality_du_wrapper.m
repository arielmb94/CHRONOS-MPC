function mpc = equality_du_wrapper(mpc)
[mpc.rp_0,mpc.rp_k] = equality_du_local(mpc.rp_0,mpc.rp_k,mpc.u,mpc.su,...
    mpc.su_col,mpc.N);
end

function [rp_0,rp_k] = equality_du_local(rp_0,rp_k,u,su,su_idx,N)
% this is the control input delay equality condition: su+ = u_prev
% [I -I][u su+]' = 0
% k = 0
rp_0(su_idx) = u(:,1)-su(:,1);
% k = 1...N
rp_k(su_idx,:) = u(:,2:N)-su(:,2:N);
end
