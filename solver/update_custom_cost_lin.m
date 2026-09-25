function mpc = update_custom_cost_lin(mpc,update_customcost_lin)
[mpc.update_customcost_lin,mpc.grad_u_Zlin_0,mpc.grad_s_Zlin,mpc.grad_su_Zlin, ...
    mpc.grad_u_Zlin,mpc.grad_s_Zlin_ter] = update_custom_cost_lin_local(update_customcost_lin, ...
    mpc.z_use_k0,mpc.z_use_s,mpc.z_use_su,mpc.z_use_u,mpc.z_use_ter, ...
    mpc.grad_u_Zlin_0,mpc.grad_s_Zlin,mpc.grad_su_Zlin,mpc.grad_u_Zlin, ...
    mpc.grad_s_Zlin_ter,mpc.qz_0,mpc.qz,mpc.qz_ter,mpc.Dz_0,mpc.Cz, ...
    mpc.Dsuz,mpc.Dz,mpc.Cz_ter,mpc.N);
end

function [update_customcost_lin,grad_u_Zlin_0,grad_s_Zlin,grad_su_Zlin, ...
    grad_u_Zlin,grad_s_Zlin_ter] = update_custom_cost_lin_local(update_customcost_lin, ...
    z_use_k0,z_use_s,z_use_su,z_use_u,z_use_ter,grad_u_Zlin_0,grad_s_Zlin, ...
    grad_su_Zlin,grad_u_Zlin,grad_s_Zlin_ter,qz_0,qz,qz_ter,Dz_0,Cz,Dsuz,Dz,Cz_ter,N)

update_customcost_lin = false;
if ~isempty(z_use_k0)
    grad_u_Zlin_0(:) = Dz_0'*qz_0;
end
if ~isempty(z_use_s) && ~isempty(z_use_su) && ~isempty(z_use_u)
    for k=1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
elseif ~isempty(z_use_s) && ~isempty(z_use_su)
    for k=1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
    end
elseif ~isempty(z_use_s) && ~isempty(z_use_u)
    for k=1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
elseif ~isempty(z_use_su) && ~isempty(z_use_u)
    for k=1:N-1
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
elseif ~isempty(z_use_s)
    for k=1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
    end
elseif ~isempty(z_use_su)
    for k=1:N-1
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
    end
elseif ~isempty(z_use_u)
    for k=1:N-1
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
end
if ~isempty(z_use_ter)
    grad_s_Zlin_ter(:) = Cz_ter'*qz_ter;
end
end
