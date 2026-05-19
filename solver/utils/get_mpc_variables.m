function mpc = get_mpc_variables(mpc,x,s_prev,u_prev,r,d,dh,dz)

%states
mpc = get_mpc_x(x,s_prev,mpc);

% control actions
mpc = get_mpc_u(x,mpc);

% differential control action
if ~isempty(mpc.du_cnstr) || ~isempty(mpc.Rdu)
    mpc = get_mpc_diff_u(u_prev,mpc);
end

if ~isempty(mpc.C)
    % system outputs
    mpc.y(:) = get_mpc_lin_out(mpc.s_all,mpc.u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,...
                        mpc.N_ctr_hor,mpc.Ny,mpc.C,mpc.D,mpc.Dd);
end

if ~isempty(r)
    % error signal
    mpc = get_error(r,mpc.y,mpc);
end

% General Linear Inequalities box constraints
if ~isempty(mpc.h_cnstr)
    % general constraints
    mpc.h(:) = get_mpc_lin_out(mpc.s_all,mpc.u,dh,mpc.nx,mpc.nu,mpc.nh,mpc.ndh,mpc.N_ctr_hor,...
        mpc.Nh,mpc.Ch,mpc.Dh,mpc.Ddh);

end

if ~isempty(mpc.Qz) || ~isempty(mpc.qz)
    % compute vector z
    mpc.z(:) = get_mpc_lin_out(mpc.s_all,mpc.u,dz,mpc.nx,mpc.nu,mpc.nz,mpc.ndz,...
        mpc.N_ctr_hor,mpc.Nz,mpc.Cz,mpc.Dz,mpc.Ddz);
end

% get slack variables 
if mpc.Nv
   mpc = get_mpc_v(x,mpc);
end

end