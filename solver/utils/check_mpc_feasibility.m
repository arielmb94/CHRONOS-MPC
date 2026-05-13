function [mpc,feas] = check_mpc_feasibility(mpc,x_ref)

feas = 1;
% State box constraints
if ~isempty(mpc.s_cnstr)

    mpc.s_cnstr = fi_box_compute(mpc.s_cnstr,mpc.s,mpc.Nx,mpc.nx);
    mpc.s_cnstr = fi_box_slack_compute(mpc,mpc.s_cnstr,mpc.Nx,mpc.nx);

    feas_box = fi_box_is_feasible(mpc.s_cnstr);
    feas_slack = fi_slack_is_feasible(mpc.s_cnstr);

    feas = feas_box && feas_slack;
end
    
% Control box constraints
if feas && ~isempty(mpc.u_cnstr)

    mpc.u_cnstr = fi_box_compute(mpc.u_cnstr,mpc.u,mpc.Nu,mpc.nu);
    feas = fi_box_is_feasible(mpc.u_cnstr);
end

% Differential Control box constraints
if feas && ~isempty(mpc.du_cnstr)
    
    mpc.du_cnstr = fi_box_compute(mpc.du_cnstr,mpc.du,mpc.Nu,mpc.nu);
    feas = fi_box_is_feasible(mpc.du_cnstr);
end

% Outputs box constraints
if feas && ~isempty(mpc.y_cnstr)

    mpc.y_cnstr = fi_box_compute(mpc.y_cnstr,mpc.y,mpc.Ny,mpc.ny);
    mpc.y_cnstr = fi_box_slack_compute(mpc,mpc.y_cnstr,mpc.Ny,mpc.ny);

    feas_box = fi_box_is_feasible(mpc.y_cnstr);
    feas_slack = fi_slack_is_feasible(mpc.y_cnstr);

    feas = feas_box && feas_slack;
end

% General Linear Inequalities box constraints
if feas && ~isempty(mpc.h_cnstr)

    mpc.h_cnstr = fi_box_compute(mpc.h_cnstr,mpc.h,mpc.Nh,mpc.nh);
    mpc.h_cnstr = fi_box_slack_compute(mpc,mpc.h_cnstr,mpc.Nh,mpc.nh);

    feas_box = fi_box_is_feasible(mpc.h_cnstr);
    feas_slack = fi_slack_is_feasible(mpc.h_cnstr);

    feas = feas_box && feas_slack;
end

% Terminal Constraint
% CODEGEN NOTE: to be commented out if there is not terminal ingredients
if feas && mpc.ter_ingredients && mpc.ter_constraint
    [mpc,feas] = get_terConst_feas(mpc,x_ref);
end

end