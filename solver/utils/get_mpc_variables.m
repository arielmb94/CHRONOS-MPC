function mpc = get_mpc_variables(mpc,has_du,tracking_cost,has_y_cnstr, ...
    has_h_cnstr,quad_custom_cost,lin_custom_cost,s_prev,u_prev)

mpc.s_ter = mpc_terminal_state_local(mpc.s_ter,mpc.s,mpc.N);
if ~isempty(has_du)
    mpc.du = mpc_diff_u_local(mpc.du,mpc.u,mpc.su,u_prev,mpc.N);
end
if ~isempty(tracking_cost) || ~isempty(has_y_cnstr)
    [mpc.y_0,mpc.y,mpc.y_ter,mpc.err_0,mpc.err,mpc.err_ter] = mpc_y_local( ...
        mpc.y_0,mpc.y,mpc.y_ter,mpc.err_0,mpc.err,mpc.err_ter,mpc.s,mpc.u,mpc.d, ...
        mpc.r_0,mpc.r,mpc.r_ter,s_prev,mpc.N,tracking_cost,mpc.y_use_k0, ...
        mpc.y_use_s,mpc.y_use_u,mpc.y_use_d,mpc.y_use_ter,mpc.C_0,mpc.D_0, ...
        mpc.Dd_0,mpc.C,mpc.D,mpc.Dd,mpc.C_ter);
end
if ~isempty(has_h_cnstr)
    mpc = mpc_h_wrapper(mpc,mpc.h_cnstr,s_prev,u_prev);
end
if ~isempty(quad_custom_cost) || ~isempty(lin_custom_cost)
    [mpc.z_0,mpc.z,mpc.z_ter] = mpc_z_local(mpc.z_0,mpc.z,mpc.z_ter,mpc.s,mpc.u, ...
        mpc.su,mpc.dz,s_prev,u_prev,mpc.N,mpc.z_use_k0,mpc.z_use_s,mpc.z_use_u, ...
        mpc.z_use_su,mpc.z_use_d,mpc.z_use_ter,mpc.Cz_0,mpc.Dz_0,mpc.Dsuz_0, ...
        mpc.Ddz_0,mpc.Cz,mpc.Dz,mpc.Dsuz,mpc.Ddz,mpc.Cz_ter);
end

end

function s_ter = mpc_terminal_state_local(s_ter,s,N)
s_ter(:) = s(:,N);
end

function du = mpc_diff_u_local(du,u,su,u_prev,N)
du(:,1) = u(:,1)-u_prev;
% delta u needs to be computed with su to match gradient definition
du(:,2:N) = u(:,2:N)-su(:,1:N-1);
end

function [y_0,y,y_ter,err_0,err,err_ter] = mpc_y_local( ...
    y_0,y,y_ter,err_0,err,err_ter,s,u,d,r_0,r,r_ter,s_prev,N, ...
    tracking_cost,y_use_k0,y_use_s,y_use_u,y_use_d,y_use_ter, ...
    C_0,D_0,Dd_0,C,D,Dd,C_ter)
if ~isempty(y_use_k0)
    y_0(:) = D_0*u(:,1);
    if ~isempty(y_use_s), y_0(:) = y_0(:) + C_0*s_prev; end
    if ~isempty(y_use_d), y_0(:) = y_0(:) + Dd_0*d(:,1); end
    if ~isempty(tracking_cost), err_0(:) = r_0 - y_0; end
end
y(:,:) = 0;
for k = 1:N-1
    if ~isempty(y_use_s), y(:,k) = y(:,k) + C(:,:,k)*s(:,k); end
    if ~isempty(y_use_u), y(:,k) = y(:,k) + D(:,:,k)*u(:,k+1); end
    if ~isempty(y_use_d), y(:,k) = y(:,k) + Dd(:,:,k)*d(:,k+1); end
end
if ~isempty(tracking_cost), err(:,:) = r-y; end
if ~isempty(y_use_ter)
    y_ter(:) = C_ter*s(:,N);
    if ~isempty(tracking_cost), err_ter(:) = r_ter - y_ter; end
end
end

function [z_0,z,z_ter] = mpc_z_local(z_0,z,z_ter,s,u,su,dz,s_prev,u_prev,N, ...
    z_use_k0,z_use_s,z_use_u,z_use_su,z_use_d,z_use_ter, ...
    Cz_0,Dz_0,Dsuz_0,Ddz_0,Cz,Dz,Dsuz,Ddz,Cz_ter)
if ~isempty(z_use_k0)
    z_0(:) = Dz_0*u(:,1);
    if ~isempty(z_use_s), z_0(:) = z_0 + Cz_0*s_prev; end
    if ~isempty(z_use_su), z_0(:) = z_0 + Dsuz_0*u_prev; end
    if ~isempty(z_use_d), z_0(:) = z_0 + Ddz_0*dz(:,1); end
end
z(:,:) = 0;
for k = 1:N-1
    if ~isempty(z_use_s), z(:,k) = z(:,k) + Cz(:,:,k)*s(:,k); end
    if ~isempty(z_use_u), z(:,k) = z(:,k) + Dz(:,:,k)*u(:,k+1); end
    if ~isempty(z_use_su), z(:,k) = z(:,k) + Dsuz(:,:,k)*su(:,k); end
    if ~isempty(z_use_d), z(:,k) = z(:,k) + Ddz(:,:,k)*dz(:,k+1); end
end
if ~isempty(z_use_ter), z_ter(:) = Cz_ter*s(:,N); end
end
