function mpc = grad_f0_scale_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter]=grad_f0_scale_local(mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter,mpc.t);
end
function [grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter] = grad_f0_scale_local(...
    grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter,t)
grad_u_f0_0 = t*grad_u_f0_0;
grad_xe_f0_k = t*grad_xe_f0_k;
grad_u_f0_k = t*grad_u_f0_k;
grad_xe_f0_ter = t*grad_xe_f0_ter;
end
