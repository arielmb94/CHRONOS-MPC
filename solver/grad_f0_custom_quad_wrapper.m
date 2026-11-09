function mpc=grad_f0_custom_quad_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter]=grad_f0_custom_quad_local(mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter,mpc.z_use_k0,mpc.z_use_ter,mpc.z_use_x,mpc.z_use_xu,mpc.z_use_u,mpc.grad_u_Z_0,mpc.grad_x_Z,mpc.grad_xu_Z,mpc.grad_u_Z,mpc.grad_x_Z_ter,mpc.z_0,mpc.z,mpc.z_ter,mpc.N,mpc.x_col,mpc.xu_col);
end
function [grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter] = grad_f0_custom_quad_local(...
                    grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter,...
                    z_use_k0,z_use_ter,z_use_x,z_use_xu,z_use_u,...
                    grad_u_Z_0,grad_x_Z,grad_xu_Z,grad_u_Z,grad_x_Z_ter,...
                    z_0,z,z_ter,N,x_idx,xu_idx)

    if ~isempty(z_use_k0)
        grad_u_f0_0 = grad_u_f0_0 + grad_u_Z_0*z_0;
    end

    for k = 1:N-1
        if ~isempty(z_use_x) && ~isempty(z_use_xu) && ~isempty(z_use_u)

            z_k = z(:,k);
            grad_xe_f0_k(x_idx,k) = grad_xe_f0_k(x_idx,k) + grad_x_Z(:,:,k)*z_k;
            grad_xe_f0_k(xu_idx,k) = grad_xe_f0_k(xu_idx,k) + grad_xu_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_x) && ~isempty(z_use_xu)

            z_k = z(:,k);
            grad_xe_f0_k(x_idx,k) = grad_xe_f0_k(x_idx,k) + grad_x_Z(:,:,k)*z_k;
            grad_xe_f0_k(xu_idx,k) = grad_xe_f0_k(xu_idx,k) + grad_xu_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_x) && ~isempty(z_use_u)

            z_k = z(:,k);
            grad_xe_f0_k(x_idx,k) = grad_xe_f0_k(x_idx,k) + grad_x_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_xu) && ~isempty(z_use_u)

            z_k = z(:,k);
            grad_xe_f0_k(xu_idx,k) = grad_xe_f0_k(xu_idx,k) + grad_xu_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_x)

            grad_xe_f0_k(x_idx,k) = grad_xe_f0_k(x_idx,k) + grad_x_Z(:,:,k)*z(:,k);
elseif ~isempty(z_use_xu)

            grad_xe_f0_k(xu_idx,k) = grad_xe_f0_k(xu_idx,k) + grad_xu_Z(:,:,k)*z(:,k);
        elseif ~isempty(z_use_u)

            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z(:,k);
        end
    end

    if ~isempty(z_use_ter)
        grad_xe_f0_ter(x_idx) = grad_xe_f0_ter(x_idx) + grad_x_Z_ter*z_ter;
    end
end
