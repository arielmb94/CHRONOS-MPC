function mpc=grad_f0_custom_lin_wrapper(mpc)
[mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,mpc.grad_xe_f0_ter]=...
    grad_f0_custom_lin_local(mpc.grad_u_f0_0,mpc.grad_xe_f0_k,mpc.grad_u_f0_k,...
    mpc.grad_xe_f0_ter,mpc.z_use_k0,mpc.z_use_ter,mpc.z_use_x,mpc.z_use_xu,mpc.z_use_u,...
    mpc.grad_u_Zlin_0,mpc.grad_x_Zlin,mpc.grad_xu_Zlin,mpc.grad_u_Zlin,mpc.grad_x_Zlin_ter,mpc.x_col,mpc.xu_col);
end
function [grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter] = grad_f0_custom_lin_local(...
                    grad_u_f0_0,grad_xe_f0_k,grad_u_f0_k,grad_xe_f0_ter,...
                    z_use_k0,z_use_ter,z_use_x,z_use_xu,z_use_u,...
                    grad_u_Z_0,grad_x_Z,grad_xu_Z,grad_u_Z,grad_x_Z_ter,...
                    x_idx,xu_idx)

if ~isempty(z_use_k0)
    grad_u_f0_0 = grad_u_f0_0 + grad_u_Z_0;
end

% k = 1:N-1
if ~isempty(z_use_x) && ~isempty(z_use_xu) && ~isempty(z_use_u)

    grad_xe_f0_k(x_idx,:) = grad_xe_f0_k(x_idx,:) + grad_x_Z;
    grad_xe_f0_k(xu_idx,:) = grad_xe_f0_k(xu_idx,:) + grad_xu_Z;
    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
elseif ~isempty(z_use_x) && ~isempty(z_use_xu)

    grad_xe_f0_k(x_idx,:) = grad_xe_f0_k(x_idx,:) + grad_x_Z;
    grad_xe_f0_k(xu_idx,:) = grad_xe_f0_k(xu_idx,:) + grad_xu_Z;
elseif ~isempty(z_use_x) && ~isempty(z_use_u)

    grad_xe_f0_k(x_idx,:) = grad_xe_f0_k(x_idx,:) + grad_x_Z;
    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
elseif ~isempty(z_use_xu) && ~isempty(z_use_u)

    grad_xe_f0_k(xu_idx,:) = grad_xe_f0_k(xu_idx,:) + grad_xu_Z;
    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
elseif ~isempty(z_use_x)

    grad_xe_f0_k(x_idx,:) = grad_xe_f0_k(x_idx,:) + grad_x_Z;
elseif ~isempty(z_use_xu)

    grad_xe_f0_k(xu_idx,:) = grad_xe_f0_k(xu_idx,:) + grad_xu_Z;
elseif ~isempty(z_use_u)

    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
end

if ~isempty(z_use_ter)
    grad_xe_f0_ter(x_idx) = grad_xe_f0_ter(x_idx) + grad_x_Z_ter;
end

end
