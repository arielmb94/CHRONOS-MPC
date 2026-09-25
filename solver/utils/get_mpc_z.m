% Computes vectors: z = Cz*s + Dz*u + Dsuz*su + Ddz*dz
function [z_0,z,z_ter] = get_mpc_z( ...
    z_0,z,z_ter,s,u,su,dz,s_prev,u_prev,N, ...
    z_use_k0,z_use_s,z_use_u,z_use_su,z_use_d,z_use_ter, ...
    Cz_0,Dz_0,Dsuz_0,Ddz_0,Cz,Dz,Dsuz,Ddz,Cz_ter)

% k=0 (only if z depends on u)
if z_use_k0
    z_0(:) = Dz_0*u(:,1);
    if z_use_s
        z_0(:) = z_0 + Cz_0*s_prev;
    end
    if z_use_su
        z_0(:) = z_0 + Dsuz_0*u_prev;
    end
    if z_use_d
        z_0(:) = z_0 + Ddz_0*dz(:,1);
    end
end

z(:,:) = 0;
for k = 1:N-1
    if z_use_s
        z(:,k) = z(:,k) + Cz(:,:,k)*s(:,k);
    end
    if z_use_u
        z(:,k) = z(:,k) + Dz(:,:,k)*u(:,k+1);
    end
    if z_use_su
        z(:,k) = z(:,k) + Dsuz(:,:,k)*su(:,k);
    end
    if z_use_d
        z(:,k) = z(:,k) + Ddz(:,:,k)*dz(:,k+1);
    end
end

if z_use_ter
    z_ter(:) = Cz_ter*s(:,N);
end

end
