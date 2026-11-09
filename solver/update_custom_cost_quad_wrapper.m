function mpc = update_custom_cost_quad_wrapper(mpc,update_customcost_quad)
[mpc.update_customcost_quad,mpc.recompute_cost_hess,mpc.grad_u_Z_0,mpc.grad_x_Z, ...
 mpc.grad_xu_Z,mpc.grad_u_Z,mpc.grad_x_Z_ter,mpc.R_Z_0,mpc.Q_Z,mpc.R_Z,mpc.Y_Z,mpc.Q_Z_ter] = ...
 update_custom_cost_quad_local(update_customcost_quad,mpc.recompute_cost_hess,mpc.grad_u_Z_0,mpc.grad_x_Z, ...
 mpc.grad_xu_Z,mpc.grad_u_Z,mpc.grad_x_Z_ter,mpc.R_Z_0,mpc.Q_Z,mpc.R_Z,mpc.Y_Z,mpc.Q_Z_ter, ...
 mpc.z_use_k0,mpc.z_use_x,mpc.z_use_xu,mpc.z_use_u,mpc.z_use_ter,mpc.Qz_0,mpc.Qz,mpc.Qz_ter, ...
 mpc.Dz_0,mpc.Cz,mpc.Dxuz,mpc.Dz,mpc.Cz_ter,mpc.x_col,mpc.xu_col,mpc.N);
end

function [update_customcost_quad,recompute_cost_hess,grad_u_Z_0,grad_x_Z,grad_xu_Z,grad_u_Z,grad_x_Z_ter, ...
 R_Z_0,Q_Z,R_Z,Y_Z,Q_Z_ter] = update_custom_cost_quad_local(update_customcost_quad,recompute_cost_hess, ...
 grad_u_Z_0,grad_x_Z,grad_xu_Z,grad_u_Z,grad_x_Z_ter,R_Z_0,Q_Z,R_Z,Y_Z,Q_Z_ter, ...
 z_use_k0,z_use_x,z_use_xu,z_use_u,z_use_ter,Qz_0,Qz,Qz_ter,Dz_0,Cz,Dxuz,Dz,Cz_ter,x_col,xu_col,N)
update_customcost_quad = false;
recompute_cost_hess = true;
if ~isempty(z_use_k0)
    grad_u_Z_0(:,:) = Dz_0'*Qz_0;
    R_Z_0(:,:) = grad_u_Z_0*Dz_0;
end
if ~isempty(z_use_x) && ~isempty(z_use_xu) && ~isempty(z_use_u)
    for k=1:N-1
        grad_x_Z(:,:,k)=Cz(:,:,k)'*Qz(:,:,k);
        grad_xu_Z(:,:,k)=Dxuz(:,:,k)'*Qz(:,:,k);
        grad_u_Z(:,:,k)=Dz(:,:,k)'*Qz(:,:,k);
        Q_Z(x_col,x_col,k)=grad_x_Z(:,:,k)*Cz(:,:,k);
        Q_Z(x_col,xu_col,k)=grad_x_Z(:,:,k)*Dxuz(:,:,k);
        Q_Z(xu_col,x_col,k)=grad_xu_Z(:,:,k)*Cz(:,:,k);
        Q_Z(xu_col,xu_col,k)=grad_xu_Z(:,:,k)*Dxuz(:,:,k);
        Y_Z(:,x_col,k)=grad_u_Z(:,:,k)*Cz(:,:,k);
        Y_Z(:,xu_col,k)=grad_u_Z(:,:,k)*Dxuz(:,:,k);
        R_Z(:,:,k)=grad_u_Z(:,:,k)*Dz(:,:,k);
    end
elseif ~isempty(z_use_x) && ~isempty(z_use_xu)
    for k=1:N-1
        grad_x_Z(:,:,k)=Cz(:,:,k)'*Qz(:,:,k);
        grad_xu_Z(:,:,k)=Dxuz(:,:,k)'*Qz(:,:,k);
        Q_Z(x_col,x_col,k)=grad_x_Z(:,:,k)*Cz(:,:,k);
        Q_Z(x_col,xu_col,k)=grad_x_Z(:,:,k)*Dxuz(:,:,k);
        Q_Z(xu_col,x_col,k)=grad_xu_Z(:,:,k)*Cz(:,:,k);
        Q_Z(xu_col,xu_col,k)=grad_xu_Z(:,:,k)*Dxuz(:,:,k);
    end
elseif ~isempty(z_use_x) && ~isempty(z_use_u)
    for k=1:N-1
        grad_x_Z(:,:,k) = Cz(:,:,k)'*Qz(:,:,k);
        grad_u_Z(:,:,k) = Dz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_x_Z(:,:,k)*Cz(:,:,k);
        Y_Z(:,:,k) = grad_u_Z(:,:,k)*Cz(:,:,k);
        R_Z(:,:,k) = grad_u_Z(:,:,k)*Dz(:,:,k);
    end
elseif ~isempty(z_use_xu) && ~isempty(z_use_u)
    for k=1:N-1
        grad_xu_Z(:,:,k) = Dxuz(:,:,k)'*Qz(:,:,k);
        grad_u_Z(:,:,k) = Dz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_xu_Z(:,:,k)*Dxuz(:,:,k);
        Y_Z(:,:,k) = grad_u_Z(:,:,k)*Dxuz(:,:,k);
        R_Z(:,:,k) = grad_u_Z(:,:,k)*Dz(:,:,k);
    end
elseif ~isempty(z_use_x)
    for k=1:N-1
        grad_x_Z(:,:,k) = Cz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_x_Z(:,:,k)*Cz(:,:,k);
    end
elseif ~isempty(z_use_xu)
    for k=1:N-1
        grad_xu_Z(:,:,k) = Dxuz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_xu_Z(:,:,k)*Dxuz(:,:,k);
    end
elseif ~isempty(z_use_u)
    for k=1:N-1
        grad_u_Z(:,:,k) = Dz(:,:,k)'*Qz(:,:,k);
        R_Z(:,:,k) = grad_u_Z(:,:,k)*Dz(:,:,k);
    end
end
if ~isempty(z_use_ter)
    grad_x_Z_ter(:,:)=Cz_ter'*Qz_ter;
    Q_Z_ter(:,:)=grad_x_Z_ter*Cz_ter;
end
end
