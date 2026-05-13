function [mpc,x0] = mpc_slack_update(mpc,x0,x_ref)

% State box constraints
if ~isempty(mpc.s_cnstr)
    [mpc.s_cnstr,v] = fi_box_slack_update(mpc,mpc.s_cnstr,mpc.s,mpc.Nx,mpc.nx);
    mpc.v(:) = v;
end

% Outputs box constraints
if ~isempty(mpc.y_cnstr)
    [mpc.y_cnstr,v] = fi_box_slack_update(mpc,mpc.y_cnstr,mpc.y,mpc.Ny,mpc.ny);
    mpc.v(:) = v;
end

% General Linear Inequalities box constraints
if ~isempty(mpc.h_cnstr)
    [mpc.h_cnstr,v] = fi_box_slack_update(mpc,mpc.h_cnstr,mpc.h,mpc.Nh,mpc.nh);
    mpc.v(:) = v;
end

% Terminal Constraint
% CODEGEN NOTE: to be commented out if there is not terminal ingredients
if mpc.ter_ingredients && mpc.ter_constraint
    mpc = get_terConst_slack_update(mpc,x_ref);
end

% update primal vector with updated slack values
x0(mpc.Nx+mpc.Nu+1:mpc.Nx+mpc.Nu+mpc.Nv) = mpc.v;

end