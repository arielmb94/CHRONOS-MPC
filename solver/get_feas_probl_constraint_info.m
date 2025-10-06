function [mpc,feas] = ...
    get_feas_probl_constraint_info(x,s_prev,u_prev,x_ref,d,dh,v,mpc)

%states
mpc = get_mpc_x(x,s_prev,mpc);

% control actions
mpc = get_mpc_u(x,mpc);

% differential control action
mpc = get_mpc_diff_u(u_prev,mpc);

% system outputs
mpc.y(:) = get_mpc_lin_out(mpc.s_all,mpc.u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N_ctr_hor,...
    mpc.Ny,mpc.C,mpc.D,mpc.Dd);

feas = 1;
% State box constraints
if ~isempty(mpc.s_cnstr)
    mpc.s_cnstr = fi_box_fun(mpc.s_cnstr,mpc.s,mpc.Nx,mpc.nx,v);
    feas = fi_box_is_feasible(mpc.s_cnstr);
end

% Terminal State box constraints
if feas && ~isempty(mpc.s_ter_cnstr)
    mpc.s_ter_cnstr = fi_box_fun(mpc.s_ter_cnstr,mpc.s_ter,mpc.nx,mpc.nx,v);
    feas = fi_box_is_feasible(mpc.s_ter_cnstr);
end
    
% Control box constraints
if feas && ~isempty(mpc.u_cnstr)
    mpc.u_cnstr = fi_box_fun(mpc.u_cnstr,mpc.u,mpc.Nu,mpc.nu,v);
    feas = fi_box_is_feasible(mpc.u_cnstr);
end

% Differential Control box constraints
if feas && ~isempty(mpc.du_cnstr)
    mpc.du_cnstr = fi_box_fun(mpc.du_cnstr,mpc.du,mpc.Nu,mpc.nu,v);
    feas = fi_box_is_feasible(mpc.du_cnstr);
end

% Outputs box constraints
if feas && ~isempty(mpc.y_cnstr)
    mpc.y_cnstr = fi_box_fun(mpc.y_cnstr,mpc.y,mpc.Ny,mpc.ny,v);
    feas = fi_box_is_feasible(mpc.y_cnstr);
end

% General Linear Inequalities box constraints
if feas && ~isempty(mpc.h_cnstr)
    % general constraints
    mpc.h(:) = get_mpc_lin_out(mpc.s_all,mpc.u,dh,mpc.nx,mpc.nu,mpc.nh,mpc.ndh,mpc.N_ctr_hor,...
        mpc.Nh,mpc.Ch,mpc.Dh,mpc.Ddh);
        
    mpc.h_cnstr = fi_box_fun(mpc.h_cnstr,mpc.h,mpc.Nh,mpc.nh,v);
    feas = fi_box_is_feasible(mpc.h_cnstr);
end

% Terminal Constraint
if feas && mpc.ter_ingredients && mpc.ter_constraint
    mpc.fi_ter_x0(:) = get_terConst_val(x_ref,mpc.s_ter,mpc.P,v);
    if mpc.fi_ter_x0>=0 
        feas = 0;
    end
end

end