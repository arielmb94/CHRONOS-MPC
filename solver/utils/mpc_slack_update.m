function [mpc,x0,feas] = mpc_slack_update(mpc,x0,x_ref)

feas = 1;
% State box constraints
if ~isempty(mpc.s_cnstr)
    [mpc.s_cnstr,v] = fi_box_slack_update(mpc,mpc.s_cnstr,mpc.s,mpc.Nx,mpc.nx);
    mpc.v(:) = v;
    feas = fi_box_is_feasible(mpc.s_cnstr);
end

% Terminal State box constraints
if feas && ~isempty(mpc.s_ter_cnstr)
    [mpc.s_ter_cnstr,v] = fi_box_slack_update(mpc,mpc.s_ter_cnstr,mpc.s_ter,mpc.nx,mpc.nx);
    mpc.v(:) = v;
    feas = fi_box_is_feasible(mpc.s_ter_cnstr);
end
    
% Control box constraints
if feas && ~isempty(mpc.u_cnstr)
    [mpc.u_cnstr,v] = fi_box_slack_update(mpc,mpc.u_cnstr,mpc.u,mpc.Nu,mpc.nu);
    mpc.v(:) = v;
    feas = fi_box_is_feasible(mpc.u_cnstr);
end

% Differential Control box constraints
if feas && ~isempty(mpc.du_cnstr)
    [mpc.du_cnstr,v] = fi_box_slack_update(mpc,mpc.du_cnstr,mpc.du,mpc.Nu,mpc.nu);
    mpc.v(:) = v;
    feas = fi_box_is_feasible(mpc.du_cnstr);
end

% Outputs box constraints
if feas && ~isempty(mpc.y_cnstr)
    [mpc.y_cnstr,v] = fi_box_slack_update(mpc,mpc.y_cnstr,mpc.y,mpc.Ny,mpc.ny);
    mpc.v(:) = v;
    feas = fi_box_is_feasible(mpc.y_cnstr);
end

% General Linear Inequalities box constraints
if feas && ~isempty(mpc.h_cnstr)
    [mpc.h_cnstr,v] = fi_box_slack_update(mpc,mpc.h_cnstr,mpc.h,mpc.Nh,mpc.nh);
    mpc.v(:) = v;
    feas = fi_box_is_feasible(mpc.h_cnstr);
end

% Terminal Constraint
% CODEGEN NOTE: to be commented out if there is not terminal ingredients
if feas && mpc.ter_ingredients && mpc.ter_constraint
    mpc = get_terConst_val(mpc,x_ref,mpc.s_ter,mpc.P,0);
end

% update primal vector with updated slack values
x0(mpc.Nx+mpc.Nu+1:mpc.Nx+mpc.Nu+mpc.Nv) = mpc.v;

end