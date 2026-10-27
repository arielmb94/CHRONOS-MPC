function mpc = f0_hess_custom_quad_wrapper(mpc)
[mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Q_f0_ter] = f0_hess_custom_quad_local(mpc.R_f0_0,mpc.Q_f0_k,mpc.R_f0_k,mpc.Y_f0_k,mpc.Q_f0_ter,mpc.z_use_k0,mpc.z_use_x,mpc.z_use_xu,mpc.z_use_u,mpc.z_use_ter,mpc.R_Z_0,mpc.Q_Z,mpc.R_Z,mpc.Y_Z,mpc.Q_Z_ter,mpc.x_col,mpc.xu_col);
end

function [R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter] = f0_hess_custom_quad_local(R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter,z_use_k0,z_use_x,z_use_xu,z_use_u,z_use_ter,R_Z_0,Q_Z,R_Z,Y_Z,Q_Z_ter,x_col,xu_col)
if ~isempty(z_use_k0)
    R_f0_0(:,:) = R_f0_0 + R_Z_0;
end
if ~isempty(z_use_x) && ~isempty(z_use_xu) && ~isempty(z_use_u)
    Q_f0_k(:,:,:) = Q_f0_k + Q_Z;
    R_f0_k(:,:,:) = R_f0_k + R_Z;
    Y_f0_k(:,:,:) = Y_f0_k + Y_Z;
elseif ~isempty(z_use_x) && ~isempty(z_use_xu)
    Q_f0_k(:,:,:) = Q_f0_k + Q_Z;
elseif ~isempty(z_use_x) && ~isempty(z_use_u)
    Q_f0_k(x_col,x_col,:) = Q_f0_k(x_col,x_col,:) + Q_Z;
    R_f0_k(:,:,:) = R_f0_k + R_Z;
    Y_f0_k(:,x_col,:) = Y_f0_k(:,x_col,:) + Y_Z;
elseif ~isempty(z_use_xu) && ~isempty(z_use_u)
    Q_f0_k(xu_col,xu_col,:) = Q_f0_k(xu_col,xu_col,:) + Q_Z;
    R_f0_k(:,:,:) = R_f0_k + R_Z;
    Y_f0_k(:,xu_col,:) = Y_f0_k(:,xu_col,:) + Y_Z;
elseif ~isempty(z_use_x)
    Q_f0_k(x_col,x_col,:) = Q_f0_k(x_col,x_col,:) + Q_Z;
elseif ~isempty(z_use_xu)
    Q_f0_k(xu_col,xu_col,:) = Q_f0_k(xu_col,xu_col,:) + Q_Z;
elseif ~isempty(z_use_u)
    R_f0_k(:,:,:) = R_f0_k + R_Z;
end
if ~isempty(z_use_ter)
    Q_f0_ter(x_col,x_col) = Q_f0_ter(x_col,x_col) + Q_Z_ter;
end
end
