function [mpc,s,s_all,s_ter,u,du,y,h,feas] = ...
    get_feas_probl_constraint_info(x,s_prev,u_prev,x_ref,d,dh,v,mpc)

%states
[s,s_all,s_ter] = get_x(x,s_prev,mpc.nx,mpc.nu,mpc.N,mpc.N_ctr_hor,mpc.Nx);
% control actions
u = get_u(x,mpc.nx,mpc.nu,mpc.N_ctr_hor,mpc.Nu);
% differential control action
du = get_diff_u(u,u_prev,mpc.nu,mpc.N_ctr_hor,mpc.Nu);
% system outputs
y = get_lin_out(s_all,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N_ctr_hor,...
    mpc.Ny,mpc.C,mpc.D,mpc.Dd);
h = 0;

feas = 1;
% State box constraints
if ~isempty(mpc.s_cnstr)
    mpc.s_cnstr = fi_box_fun(mpc.s_cnstr,s,mpc.Nx,mpc.nx,v);
    feas = fi_box_is_feasible(mpc.s_cnstr);
end

% Terminal State box constraints
if feas && ~isempty(mpc.s_ter_cnstr)
    mpc.s_ter_cnstr = fi_box_fun(mpc.s_ter_cnstr,s_ter,mpc.nx,mpc.nx,v);
    feas = fi_box_is_feasible(mpc.s_ter_cnstr);
end
    
% Control box constraints
if feas && ~isempty(mpc.u_cnstr)
    mpc.u_cnstr = fi_box_fun(mpc.u_cnstr,u,mpc.Nu,mpc.nu,v);
    feas = fi_box_is_feasible(mpc.u_cnstr);
end

% Differential Control box constraints
if feas && ~isempty(mpc.du_cnstr)
    mpc.du_cnstr = fi_box_fun(mpc.du_cnstr,du,mpc.Nu,mpc.nu,v);
    feas = fi_box_is_feasible(mpc.du_cnstr);
end

% Outputs box constraints
if feas && ~isempty(mpc.y_cnstr)
    mpc.y_cnstr = fi_box_fun(mpc.y_cnstr,y,mpc.Ny,mpc.ny,v);
    feas = fi_box_is_feasible(mpc.y_cnstr);
end

% General Linear Inequalities box constraints
if feas && ~isempty(mpc.h_cnstr)
    % general constraints
    h = get_lin_out(s_all,u,dh,mpc.nx,mpc.nu,mpc.nh,mpc.ndh,mpc.N_ctr_hor,...
        mpc.Nh,mpc.Ch,mpc.Dh,mpc.Ddh);

    mpc.h_cnstr = fi_box_fun(mpc.h_cnstr,h,mpc.Nh,mpc.nh,v);
    feas = fi_box_is_feasible(mpc.h_cnstr);
end

% Terminal Constraint
if feas && mpc.ter_ingredients && mpc.ter_constraint
    mpc.fi_ter_x0 = get_terConst_val(x_ref,s_ter,mpc.P,v);
    if mpc.fi_ter_x0>=0 
        feas = 0;
    end
end

end