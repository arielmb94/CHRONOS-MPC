function [mpc,feas] = check_mpc_feasibility(mpc,x_ref)

feas = 1;
% State box constraints
if ~isempty(mpc.s_cnstr)
    mpc.s_cnstr = fi_box_feas_fun(mpc,mpc.s_cnstr,mpc.s,mpc.Nx,mpc.nx);
    feas = fi_box_is_feasible(mpc.s_cnstr);
end

% Terminal State box constraints
if feas && ~isempty(mpc.s_ter_cnstr)
    mpc.s_ter_cnstr = fi_box_feas_fun(mpc,mpc.s_ter_cnstr,mpc.s_ter,mpc.nx,mpc.nx);
    feas = fi_box_is_feasible(mpc.s_ter_cnstr);
end
    
% Control box constraints
if feas && ~isempty(mpc.u_cnstr)
    mpc.u_cnstr = fi_box_feas_fun(mpc,mpc.u_cnstr,mpc.u,mpc.Nu,mpc.nu);
    feas = fi_box_is_feasible(mpc.u_cnstr);
end

% Differential Control box constraints
if feas && ~isempty(mpc.du_cnstr)
    mpc.du_cnstr = fi_box_feas_fun(mpc,mpc.du_cnstr,mpc.du,mpc.Nu,mpc.nu);
    feas = fi_box_is_feasible(mpc.du_cnstr);
end

% Outputs box constraints
if feas && ~isempty(mpc.y_cnstr)
    mpc.y_cnstr = fi_box_feas_fun(mpc,mpc.y_cnstr,mpc.y,mpc.Ny,mpc.ny);
    feas = fi_box_is_feasible(mpc.y_cnstr);
end

% General Linear Inequalities box constraints
if feas && ~isempty(mpc.h_cnstr)
    mpc.h_cnstr = fi_box_feas_fun(mpc,mpc.h_cnstr,mpc.h,mpc.Nh,mpc.nh);
    feas = fi_box_is_feasible(mpc.h_cnstr);
end

% Terminal Constraint
% CODEGEN NOTE: to be commented out if there is not terminal ingredients
if feas && mpc.ter_ingredients && mpc.ter_constraint
    [mpc,feas] = get_terConst_feas(mpc,x_ref);
end

end