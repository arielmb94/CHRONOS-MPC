function mpc = f0_hess_control_wrapper(mpc)
[mpc.R_f0_0,mpc.R_f0_k] = f0_hess_control_local(mpc.R_f0_0,mpc.R_f0_k,mpc.Ru,mpc.N);
end

function [R_f0_0,R_f0_k] = f0_hess_control_local(R_f0_0,R_f0_k,Ru,N)
R_f0_0(:,:) = R_f0_0 + Ru(:,:,1);
for k = 1:N-1
    R_f0_k(:,:,k) = R_f0_k(:,:,k) + Ru(:,:,k+1);
end
end
