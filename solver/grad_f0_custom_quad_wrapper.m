function mpc=grad_f0_custom_quad_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter]=grad_f0_custom_quad_local(mpc.grad_u_f0_0,mpc.grad_se_f0_k,mpc.grad_u_f0_k,mpc.grad_se_f0_ter,mpc.z_use_k0,mpc.z_use_ter,mpc.z_use_s,mpc.z_use_su,mpc.z_use_u,mpc.grad_u_Z_0,mpc.grad_s_Z,mpc.grad_su_Z,mpc.grad_u_Z,mpc.grad_s_Z_ter,mpc.z_0,mpc.z,mpc.z_ter,mpc.N,mpc.s_col,mpc.su_col);
end
function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = grad_f0_custom_quad_local(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    z_use_k0,z_use_ter,z_use_s,z_use_su,z_use_u,...
                    grad_u_Z_0,grad_s_Z,grad_su_Z,grad_u_Z,grad_s_Z_ter,...
                    z_0,z,z_ter,N,s_idx,su_idx)

    if ~isempty(z_use_k0)
        grad_u_f0_0 = grad_u_f0_0 + grad_u_Z_0*z_0;
    end

    for k = 1:N-1
        if ~isempty(z_use_s) && ~isempty(z_use_su) && ~isempty(z_use_u)

            z_k = z(:,k);
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z_k;
            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_s) && ~isempty(z_use_su)

            z_k = z(:,k);
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z_k;
            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_s) && ~isempty(z_use_u)

            z_k = z(:,k);
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_su) && ~isempty(z_use_u)

            z_k = z(:,k);
            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif ~isempty(z_use_s)

            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z(:,k);
elseif ~isempty(z_use_su)

            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z(:,k);
        elseif ~isempty(z_use_u)

            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z(:,k);
        end
    end

    if ~isempty(z_use_ter)
        grad_se_f0_ter(s_idx) = grad_se_f0_ter(s_idx) + grad_s_Z_ter*z_ter;
    end
end
