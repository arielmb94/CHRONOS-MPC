function mpc = f0_hess_control_rate_wrapper(mpc)
[mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k] = f0_hess_control_rate_local(mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Rdu,mpc.su_col,mpc.N);
end

function [R_f0_0,Q_f0_k,R_f0_k,Y_f0_k] = f0_hess_control_rate_local(R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Rdu,su_col,N)
R_f0_0(:,:) = R_f0_0 + Rdu(:,:,1);
for k = 1:N-1
    Q_f0_k(su_col,su_col,k) = Q_f0_k(su_col,su_col,k) + Rdu(:,:,k+1);
    R_f0_k(:,:,k) = R_f0_k(:,:,k) + Rdu(:,:,k+1);
    Y_f0_k(:,su_col,k) = Y_f0_k(:,su_col,k) - Rdu(:,:,k+1);
end
end
