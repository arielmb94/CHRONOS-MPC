function mpc = equality_du_wrapper(mpc)
[mpc.rp_0,mpc.rp_k] = equality_du_local(mpc.rp_0,mpc.rp_k,mpc.u,mpc.xu,...
    mpc.xu_col,mpc.N);
end

function [rp_0,rp_k] = equality_du_local(rp_0,rp_k,u,xu,xu_idx,N)
% this is the control input delay equality condition: xu+ = u_prev
% [I -I][u xu+]' = 0
% k = 0
rp_0(xu_idx) = u(:,1)-xu(:,1);
% k = 1...N
rp_k(xu_idx,:) = u(:,2:N)-xu(:,2:N);
end
