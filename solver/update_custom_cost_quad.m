function mpc = update_custom_cost_quad(mpc,update_customcost_quad)

[mpc.update_customcost_quad,mpc.recompute_cost_hess,mpc.grad_u_Z_0,mpc.grad_s_Z, ...
 mpc.grad_su_Z,mpc.grad_u_Z,mpc.grad_s_Z_ter,mpc.R_Z_0,mpc.Q_Z,mpc.R_Z,mpc.Y_Z,mpc.Q_Z_ter] = ...
 update_custom_cost_quad_local(update_customcost_quad,mpc.recompute_cost_hess,mpc.grad_u_Z_0,mpc.grad_s_Z, ...
 mpc.grad_su_Z,mpc.grad_u_Z,mpc.grad_s_Z_ter,mpc.R_Z_0,mpc.Q_Z,mpc.R_Z,mpc.Y_Z,mpc.Q_Z_ter, ...
 mpc.z_use_k0,mpc.z_use_s,mpc.z_use_su,mpc.z_use_u,mpc.z_use_ter,mpc.Qz_0,mpc.Qz,mpc.Qz_ter, ...
 mpc.Dz_0,mpc.Cz,mpc.Dsuz,mpc.Dz,mpc.Cz_ter,mpc.s_col,mpc.su_col,mpc.N);
end

function [update_customcost_quad,recompute_cost_hess,grad_u_Z_0,grad_s_Z,grad_su_Z,grad_u_Z,grad_s_Z_ter, ...
 R_Z_0,Q_Z,R_Z,Y_Z,Q_Z_ter] = update_custom_cost_quad_local(update_customcost_quad,recompute_cost_hess, ...
 grad_u_Z_0,grad_s_Z,grad_su_Z,grad_u_Z,grad_s_Z_ter,R_Z_0,Q_Z,R_Z,Y_Z,Q_Z_ter, ...
 z_use_k0,z_use_s,z_use_su,z_use_u,z_use_ter,Qz_0,Qz,Qz_ter,Dz_0,Cz,Dsuz,Dz,Cz_ter,s_col,su_col,N)

update_customcost_quad = false;
recompute_cost_hess = true;
if ~isempty(z_use_k0)
    grad_u_Z_0(:,:) = Dz_0'*Qz_0;
    R_Z_0(:,:) = grad_u_Z_0*Dz_0;
end
if ~isempty(z_use_s) && ~isempty(z_use_su) && ~isempty(z_use_u)
    for k=1:N-1
        grad_s_Z(:,:,k)=Cz(:,:,k)'*Qz(:,:,k);
        grad_su_Z(:,:,k)=Dsuz(:,:,k)'*Qz(:,:,k);
        grad_u_Z(:,:,k)=Dz(:,:,k)'*Qz(:,:,k);
        Q_Z(s_col,s_col,k)=grad_s_Z(:,:,k)*Cz(:,:,k);
        Q_Z(s_col,su_col,k)=grad_s_Z(:,:,k)*Dsuz(:,:,k);
        Q_Z(su_col,s_col,k)=grad_su_Z(:,:,k)*Cz(:,:,k);
        Q_Z(su_col,su_col,k)=grad_su_Z(:,:,k)*Dsuz(:,:,k);
        Y_Z(:,s_col,k)=grad_u_Z(:,:,k)*Cz(:,:,k);
        Y_Z(:,su_col,k)=grad_u_Z(:,:,k)*Dsuz(:,:,k);
        R_Z(:,:,k)=grad_u_Z(:,:,k)*Dz(:,:,k);
    end
elseif ~isempty(z_use_s) && ~isempty(z_use_su)
    for k=1:N-1
        grad_s_Z(:,:,k)=Cz(:,:,k)'*Qz(:,:,k);
        grad_su_Z(:,:,k)=Dsuz(:,:,k)'*Qz(:,:,k);
        Q_Z(s_col,s_col,k)=grad_s_Z(:,:,k)*Cz(:,:,k);
        Q_Z(s_col,su_col,k)=grad_s_Z(:,:,k)*Dsuz(:,:,k);
        Q_Z(su_col,s_col,k)=grad_su_Z(:,:,k)*Cz(:,:,k);
        Q_Z(su_col,su_col,k)=grad_su_Z(:,:,k)*Dsuz(:,:,k);
    end
elseif ~isempty(z_use_s) && ~isempty(z_use_u)
    for k=1:N-1
        grad_s_Z(:,:,k) = Cz(:,:,k)'*Qz(:,:,k);
        grad_u_Z(:,:,k) = Dz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_s_Z(:,:,k)*Cz(:,:,k);
        Y_Z(:,:,k) = grad_u_Z(:,:,k)*Cz(:,:,k);
        R_Z(:,:,k) = grad_u_Z(:,:,k)*Dz(:,:,k);
    end
elseif ~isempty(z_use_su) && ~isempty(z_use_u)
    for k=1:N-1
        grad_su_Z(:,:,k) = Dsuz(:,:,k)'*Qz(:,:,k);
        grad_u_Z(:,:,k) = Dz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_su_Z(:,:,k)*Dsuz(:,:,k);
        Y_Z(:,:,k) = grad_u_Z(:,:,k)*Dsuz(:,:,k);
        R_Z(:,:,k) = grad_u_Z(:,:,k)*Dz(:,:,k);
    end
elseif ~isempty(z_use_s)
    for k=1:N-1
        grad_s_Z(:,:,k) = Cz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_s_Z(:,:,k)*Cz(:,:,k);
    end
elseif ~isempty(z_use_su)
    for k=1:N-1
        grad_su_Z(:,:,k) = Dsuz(:,:,k)'*Qz(:,:,k);
        Q_Z(:,:,k) = grad_su_Z(:,:,k)*Dsuz(:,:,k);
    end
elseif ~isempty(z_use_u)
    for k=1:N-1
        grad_u_Z(:,:,k) = Dz(:,:,k)'*Qz(:,:,k);
        R_Z(:,:,k) = grad_u_Z(:,:,k)*Dz(:,:,k);
    end
end
if ~isempty(z_use_ter)
    grad_s_Z_ter(:,:)=Cz_ter'*Qz_ter;
    Q_Z_ter(:,:)=grad_s_Z_ter*Cz_ter;
end
end

