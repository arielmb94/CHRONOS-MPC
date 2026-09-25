function mpc=grad_f0_tracking_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter]=grad_f0_tracking_local(mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter,mpc.y_use_k0,mpc.y_use_ter,mpc.y_use_s,mpc.y_use_u,mpc.grad_u_E_0,mpc.grad_s_E,mpc.grad_u_E,mpc.grad_s_E_ter,mpc.err_0,mpc.err,mpc.err_ter,mpc.N,mpc.s_col);
end
function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = grad_f0_tracking_local(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    y_use_k0,y_use_ter,y_use_s,y_use_u,...
                    grad_u_E_0,grad_s_E_k,grad_u_E_k,grad_s_E_ter,...
                    err_0,err,err_ter,N,s_idx)

    if ~isempty(y_use_k0)
        grad_u_f0_0 = grad_u_f0_0 - grad_u_E_0*err_0;
    end

    for k = 1:N-1
        if ~isempty(y_use_s) && ~isempty(y_use_u)
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) -...
                                        grad_s_E_k(:,:,k)*err(:,k);
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) -...
                                        grad_u_E_k(:,:,k)*err(:,k);
        elseif ~isempty(y_use_s)
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) -...
                                        grad_s_E_k(:,:,k)*err(:,k);
        elseif ~isempty(y_use_u)
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) -...
                                        grad_u_E_k(:,:,k)*err(:,k);
        end
    end

    if ~isempty(y_use_ter)
        grad_se_f0_ter(s_idx) = grad_se_f0_ter(s_idx) - grad_s_E_ter*err_ter;
    end
end
