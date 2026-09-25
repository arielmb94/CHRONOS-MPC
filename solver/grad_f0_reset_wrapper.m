function mpc = grad_f0_reset_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter]=...
    grad_f0_reset_local(mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter);
end

function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = grad_f0_reset_local(...
    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter)
grad_u_f0_0(:) = 0;
grad_se_f0_k(:,:) = 0;
grad_u_f0_k(:,:) = 0;
grad_se_f0_ter(:) = 0;
end
