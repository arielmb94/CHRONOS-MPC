function mpc=grad_f0_control_rate_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k]=grad_f0_control_rate_local(mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.Rdu,mpc.du,mpc.xu_col,mpc.N);
end
function [grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k] = grad_f0_control_rate_local(grad_u_f0_0,...
                                    grad_xe_f0_k,grad_u_f0_k,Rdu,du,xu_idx,N)
    % k = 0
    grad_u_f0_0 = grad_u_f0_0 + Rdu(:,:,1)*du(:,1);

    for k = 1:N-1
        grad_du_k = Rdu(:,:,k+1)*du(:,k+1);

        grad_xe_f0_k(xu_idx,k) = grad_xe_f0_k(xu_idx,k) - grad_du_k;
        grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_du_k;
    end
end
