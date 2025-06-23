function [s,s_all,s_ter,u,du,y,err,h,...
    fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
    fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
    fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
    fi_h_min_x0,fi_h_max_x0,feas] = ...
    get_state_constraint_info(x,s_prev,u_prev,r,x_ref,d,dh,mpc)

%states
[s,s_all,s_ter] = get_x(x,s_prev,mpc.nx,mpc.nu,mpc.N,mpc.N_ctr_hor,mpc.Nx);
% control actions
u = get_u(x,mpc.nx,mpc.nu,mpc.N_ctr_hor,mpc.Nu);
% differential control action
du = get_diff_u(u,u_prev,mpc.nu,mpc.N_ctr_hor,mpc.Nu);
% system outputs
y = get_lin_out(s_all,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N_ctr_hor,...
    mpc.Ny,mpc.C,mpc.D,mpc.Dd);
% error signal
err = get_error(r,y);

feas = 1;
% State box constraints
if ~isempty(mpc.x_min) || ~isempty(mpc.x_max)
    [fi_s_min_x0,fi_s_max_x0] = fi_box_fun(s,mpc.x_min,mpc.x_max,mpc.Nx,mpc.nx,0);
    if any(fi_s_min_x0>0) || any(fi_s_max_x0>0)
        feas = 0;
    end
else
    fi_s_min_x0 = []; 
    fi_s_max_x0 = [];
end

% Terminal State box constraints
if feas && (~isempty(mpc.x_ter_min) || ~isempty(mpc.x_ter_max))
    [fi_s_ter_min_x0,fi_s_ter_max_x0] = fi_box_fun(s_ter,mpc.x_ter_min,mpc.x_ter_max,mpc.nx,mpc.nx,0);
    if any(fi_s_ter_min_x0>0) || any(fi_s_ter_max_x0>0)
        feas = 0;
    end
else
    fi_s_ter_min_x0 = [];
    fi_s_ter_max_x0 = [];
end

% Control box constraints
if feas && (~isempty(mpc.u_min) || ~isempty(mpc.u_max))
    [fi_u_min_x0,fi_u_max_x0] = fi_box_fun(u,mpc.u_min,mpc.u_max,mpc.Nu,mpc.nu,0);
    if any(fi_u_min_x0>0) || any(fi_u_max_x0>0)
        feas = 0;
    end
else
    fi_u_min_x0 = [];
    fi_u_max_x0 = [];
end

% Differential Control box constraints
if feas && (~isempty(mpc.du_min) || ~isempty(mpc.du_max))
    [fi_du_min_x0,fi_du_max_x0] = fi_box_fun(du,mpc.du_min,mpc.du_max,mpc.Nu,mpc.nu,0);
    if any(fi_du_min_x0>0) || any(fi_du_max_x0>0)
        feas = 0;
    end
else
    fi_du_min_x0 = [];
    fi_du_max_x0 = [];
end

% Outputs box constraints
if feas && (~isempty(mpc.y_min) || ~isempty(mpc.y_max))
    [fi_y_min_x0,fi_y_max_x0] = fi_box_fun(y,mpc.y_min,mpc.y_max,mpc.Ny,mpc.ny,0);
    if any(fi_y_min_x0>0) || any(fi_y_max_x0>0)
        feas = 0;
    end
else
    fi_y_min_x0 = [];
    fi_y_max_x0 = [];
end

% Terminal Constraint
if feas && mpc.ter_ingredients && mpc.ter_constraint
    fi_ter_x0 = get_terConst_val(x_ref,s_ter,mpc.P,0);
    if fi_ter_x0>=0 
        feas = 0;
    end
else
    fi_ter_x0 = [];
end

% General Linear Inequalities box constraints
if feas && (~isempty(mpc.h_min) || ~isempty(mpc.h_max))
    % general constraints
    h = get_lin_out(s_all,u,dh,mpc.nx,mpc.nu,mpc.nh,mpc.ndh,mpc.N_ctr_hor,...
        mpc.Nh,mpc.Ch,mpc.Dh,mpc.Ddh);

    [fi_h_min_x0,fi_h_max_x0] = fi_box_fun(h,mpc.h_min,mpc.h_max,mpc.Nh,mpc.nh,0);
    if any(fi_h_min_x0>0) || any(fi_h_max_x0>0)
        feas = 0;
    end
else
    fi_h_min_x0 = [];
    fi_h_max_x0 = [];
    h = [];
end

end