function mpc=grad_f0_control_wrapper(mpc,quad_control_cost,lin_control_cost)
[mpc.grad_u_f0_0,mpc.grad_u_f0_k]=grad_f0_control_local(quad_control_cost,lin_control_cost,mpc.grad_u_f0_0,mpc.grad_u_f0_k,mpc.Ru,mpc.ru,mpc.u,mpc.N);
end

function [grad_u_f0_0,grad_u_f0_k] = grad_f0_control_local(quad_control_cost,lin_control_cost,...
                                    grad_u_f0_0,grad_u_f0_k,Ru,ru,u,N)
    if ~isempty(quad_control_cost)
        grad_u_f0_0 = grad_u_f0_0 + Ru(:,:,1)*u(:,1);
    
        for k = 1:N-1
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + Ru(:,:,k+1)*u(:,k+1);
        end
    end
    if ~isempty(lin_control_cost)
        grad_u_f0_0 = grad_u_f0_0 + ru(:,1);
    
        grad_u_f0_k = grad_u_f0_k + ru(:,2:N);
    end
end
