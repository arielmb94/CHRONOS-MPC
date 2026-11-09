function mpc=grad_f0_tracking_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter]=grad_f0_tracking_local(mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter,mpc.y_use_k0,mpc.y_use_ter,mpc.y_use_x,mpc.y_use_u,mpc.grad_u_E_0,mpc.grad_x_E,mpc.grad_u_E,mpc.grad_x_E_ter,mpc.err_0,mpc.err,mpc.err_ter,mpc.N,mpc.x_col);
end
function [grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter] = grad_f0_tracking_local(...
                    grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter,...
                    y_use_k0,y_use_ter,y_use_x,y_use_u,...
                    grad_u_E_0,grad_x_E_k,grad_u_E_k,grad_x_E_ter,...
                    err_0,err,err_ter,N,x_idx)

    if ~isempty(y_use_k0)
        grad_u_f0_0 = grad_u_f0_0 - grad_u_E_0*err_0;
    end

    for k = 1:N-1
        if ~isempty(y_use_x) && ~isempty(y_use_u)
            grad_xe_f0_k(x_idx,k) = grad_xe_f0_k(x_idx,k) -...
                                        grad_x_E_k(:,:,k)*err(:,k);
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) -...
                                        grad_u_E_k(:,:,k)*err(:,k);
        elseif ~isempty(y_use_x)
            grad_xe_f0_k(x_idx,k) = grad_xe_f0_k(x_idx,k) -...
                                        grad_x_E_k(:,:,k)*err(:,k);
        elseif ~isempty(y_use_u)
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) -...
                                        grad_u_E_k(:,:,k)*err(:,k);
        end
    end

    if ~isempty(y_use_ter)
        grad_xe_f0_ter(x_idx) = grad_xe_f0_ter(x_idx) - grad_x_E_ter*err_ter;
    end
end
