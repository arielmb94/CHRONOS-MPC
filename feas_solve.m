function [x_mpc,iter,mpc] = feas_solve(x0,mpc,s_prev,u_prev,d,x_ref,di)
arguments
x0
mpc
s_prev
u_prev
d = [];
x_ref = [];
di = [];
end

    % Set false at start
    mpc.unfeasible = 0;

    x_mpc = x0;
    v_feas = mpc.v0_feas;
    % extend optimization vector with feasibility slack variable
    x = [x_mpc;v_feas];

    % number of variables of original mpc problem
    n = length(x);
    
    % Generate constraint's gradient from original mpc
    feas_slv = feas_init(mpc);

    % Adapt equality constraints to feasibility solver size 
    % number of equality constraints
    n_eq = size(mpc.Aeq,1); 
    Aeq_feas = [mpc.Aeq zeros(n_eq,1)];

    % get mpc variables from optimization vector x and constraint
    % information
    % increase v_feas if needed   
    feas = 0;
    while ~feas

        [s,s_all,s_ter,u,du,y,yi,...
            fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
            fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
            fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
            fi_yi_min_x0,fi_yi_max_x0,feas] = ...
            get_feas_probl_constraint_info(x_mpc,s_prev,u_prev,x_ref,d,di,...
            v_feas,mpc);

        if ~feas
            v_feas = v_feas * mpc.feas_lambda;
            x(end) = v_feas;
        end
    end

    opts.SYM = true;
    iter = 0;
    while v_feas > 0 && iter <= mpc.max_feas_iter

        % Compute gradient:

        % 1. Compute gradient of box inequalities at x0:
        % init inequalities gradient vector
        grad_fi_Ind = zeros(n,1);

        % state inequalities
        if ~isempty(mpc.x_min)
            grad_s_min_Ind_x0 = grad_box_Ind(fi_s_min_x0,feas_slv.gradXmin);

            grad_fi_Ind = grad_fi_Ind + grad_s_min_Ind_x0;
        end

        if ~isempty(mpc.x_max)
            grad_x_max_Ind_x0 = grad_box_Ind(fi_s_max_x0,feas_slv.gradXmax);

            grad_fi_Ind = grad_fi_Ind + grad_x_max_Ind_x0;
        end

        % terminal state inequalities
        if ~isempty(mpc.x_ter_min)
            grad_s_ter_min_Ind_x0 = grad_box_Ind(fi_s_ter_min_x0,feas_slv.gradXtermin);

            grad_fi_Ind = grad_fi_Ind + grad_s_ter_min_Ind_x0;
        end

        if ~isempty(mpc.x_ter_max)
            grad_x_ter_max_Ind_x0 = grad_box_Ind(fi_s_ter_max_x0,feas_slv.gradXtermax);

            grad_fi_Ind = grad_fi_Ind + grad_x_ter_max_Ind_x0;
        end

        % control inequalities
        if ~isempty(mpc.u_min)
            grad_u_min_Ind_x0 = grad_box_Ind(fi_u_min_x0,feas_slv.gradUmin);

            grad_fi_Ind = grad_fi_Ind + grad_u_min_Ind_x0;
        end

        if ~isempty(mpc.u_max)
            grad_u_max_Ind_x0 = grad_box_Ind(fi_u_max_x0,feas_slv.gradUmax);

            grad_fi_Ind = grad_fi_Ind + grad_u_max_Ind_x0;
        end

        % control differential inequalities
        if ~isempty(mpc.du_min)
            grad_du_min_Ind_x0 = grad_box_Ind(fi_du_min_x0,feas_slv.gradDeltaUmin);

            grad_fi_Ind = grad_fi_Ind + grad_du_min_Ind_x0;
        end

        if ~isempty(mpc.du_max)
            grad_du_max_Ind_x0 = grad_box_Ind(fi_du_max_x0,feas_slv.gradDeltaUmax);

            grad_fi_Ind = grad_fi_Ind + grad_du_max_Ind_x0;
        end

        % output inequalities
        if ~isempty(mpc.y_min)
            grad_y_min_Ind_x0 = grad_box_Ind(fi_y_min_x0,feas_slv.gradYmin);

            grad_fi_Ind = grad_fi_Ind + grad_y_min_Ind_x0;
        end

        if ~isempty(mpc.y_max)
            grad_y_max_Ind_x0 = grad_box_Ind(fi_y_max_x0,feas_slv.gradYmax);

            grad_fi_Ind = grad_fi_Ind + grad_y_max_Ind_x0;
        end

        % General Linear inequalities
        if ~isempty(mpc.yi_min)
            grad_yi_min_Ind_x0 = grad_box_Ind(fi_yi_min_x0,feas_slv.gradYimin);

            grad_fi_Ind = grad_fi_Ind + grad_yi_min_Ind_x0;
        end

        if ~isempty(mpc.yi_max)
            grad_yi_max_Ind_x0 = grad_box_Ind(fi_yi_max_x0,feas_slv.gradYimax);

            grad_fi_Ind = grad_fi_Ind + grad_yi_max_Ind_x0;
        end

        % 2. If enabled, compute terminal constraint gradients 
        if mpc.ter_ingredients && mpc.ter_constraint

            [grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
                ter_set_feas_Ind_fun(x_ref,s_ter,fi_ter_x0,mpc.P,...
                mpc.nx,n);

                grad_fi_Ind = grad_fi_Ind + grad_ter_Ind_x0; 
        end


        % Compute Hessian Matrix

        % 1. Compute Hessian of box inequalities at x0:

        % init inequalities hessian vector
        hess_fi_Ind = zeros(n);

        % state inequalities
        if ~isempty(mpc.x_min)
            hess_s_min_Ind_x0 = hess_linear_Ind(fi_s_min_x0,feas_slv.hessXmin);

            hess_fi_Ind = hess_fi_Ind + hess_s_min_Ind_x0;
        end

        if ~isempty(mpc.x_max)
            hess_s_max_Ind_x0 = hess_linear_Ind(fi_s_max_x0,feas_slv.hessXmax);

            hess_fi_Ind = hess_fi_Ind + hess_s_max_Ind_x0;
        end

        % terminal state inequalities
        if ~isempty(mpc.x_ter_min)
            hess_s_ter_min_Ind_x0 = hess_linear_Ind(fi_s_ter_min_x0,feas_slv.hessXtermin);

            hess_fi_Ind = hess_fi_Ind + hess_s_ter_min_Ind_x0;
        end

        if ~isempty(mpc.x_ter_max)
            hess_s_ter_max_Ind_x0 = hess_linear_Ind(fi_s_ter_max_x0,feas_slv.hessXtermax);

            hess_fi_Ind = hess_fi_Ind + hess_s_ter_max_Ind_x0;
        end

        % control inequalities
        if ~isempty(mpc.u_min)
            hess_u_min_Ind_x0 = hess_linear_Ind(fi_u_min_x0,feas_slv.hessUmin);

            hess_fi_Ind = hess_fi_Ind + hess_u_min_Ind_x0;
        end

        if ~isempty(mpc.u_max)
            hess_u_max_Ind_x0 = hess_linear_Ind(fi_u_max_x0,feas_slv.hessUmax);

            hess_fi_Ind = hess_fi_Ind + hess_u_max_Ind_x0;
        end

        % control differential inequalities
        if ~isempty(mpc.du_min)
            hess_du_min_Ind_x0 = hess_linear_Ind(fi_du_min_x0,feas_slv.hessDeltaUmin);

            hess_fi_Ind = hess_fi_Ind + hess_du_min_Ind_x0;
        end

        if ~isempty(mpc.du_max)
            hess_du_max_Ind_x0 = hess_linear_Ind(fi_du_max_x0,feas_slv.hessDeltaUmax);

            hess_fi_Ind = hess_fi_Ind + hess_du_max_Ind_x0;
        end

        % output inequalities
        if ~isempty(mpc.y_min)
            hess_y_min_Ind_x0 = hess_linear_Ind(fi_y_min_x0,feas_slv.hessYmin);

            hess_fi_Ind = hess_fi_Ind + hess_y_min_Ind_x0;
        end

        if ~isempty(mpc.y_max)
            hess_y_max_Ind_x0 = hess_linear_Ind(fi_y_max_x0,feas_slv.hessYmax);

            hess_fi_Ind = hess_fi_Ind + hess_y_max_Ind_x0;
        end

        % General Linear inequalities
        if ~isempty(mpc.yi_min)
            hess_yi_min_Ind_x0 = hess_linear_Ind(fi_yi_min_x0,feas_slv.hessYimin);

            hess_fi_Ind = hess_fi_Ind + hess_yi_min_Ind_x0;
        end

        if ~isempty(mpc.yi_max)
            hess_yi_max_Ind_x0 = hess_linear_Ind(fi_yi_max_x0,feas_slv.hessYimax);

            hess_fi_Ind = hess_fi_Ind + hess_yi_max_Ind_x0;
        end

        % 2. If enabled, add terminal constraint hessian term
        if mpc.ter_ingredients && mpc.ter_constraint
            hess_fi_Ind = hess_fi_Ind + hess_ter_Ind_x0;
        end


        % Manage and adapt gradients and Hessian
        
        % 1. Compute gradient of f0
        [grad_f0,hess_f0] = grad_hess_f0_feas(mpc,x_mpc,x0,n-1);

        % 2. Compute gradient at x0 : grad(J) = t*grad(f0)+grad(Phi)
        grad_J_x0 = mpc.t_feas*grad_f0 + grad_fi_Ind;

        % 3. Compute Hessian of f(x0,t):
        if mpc.warm_starting

            hess_J_x0 = hess_fi_Ind;
        else  

            hess_J_x0 = mpc.t_feas*hess_f0 + hess_fi_Ind;
        end

        % solve KKT system

        KKT = [hess_J_x0 Aeq_feas';Aeq_feas zeros(n_eq)];

        delta_x = - linsolve(KKT,[grad_J_x0;Aeq_feas*x-mpc.beq],opts);
        %delta_x = - linsolve(KKT,[grad_J_x0;zeros(n_eq,1)],opts);
        %delta_x = - KKT\[grad_J_x0;Aeq_feas*x-mpc.beq];
        delta_x_prim = delta_x(1:n);

        % Feasibility line search
        l = 1;
        xhat = x+l*delta_x_prim;
        x_mpc = xhat(1:n-1);
        v_feas = xhat(n);

        [s,s_all,s_ter,u,du,y,yi,...
         fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
         fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
         fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
         fi_yi_min_x0,fi_yi_max_x0,feas] = ...
         get_feas_probl_constraint_info(x_mpc,s_prev,u_prev,x_ref,d,di,v_feas,mpc);

        if feas
            x = xhat;
        else
            while ~feas
                
                l = l*mpc.Beta;

                xhat = x+l*delta_x_prim;
                x_mpc = xhat(1:n-1);
                v_feas = xhat(n);

                [s,s_all,s_ter,u,du,y,yi,...
                 fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
                 fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
                 fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
                 fi_yi_min_x0,fi_yi_max_x0,feas] = ...
                 get_feas_probl_constraint_info(x_mpc,s_prev,u_prev,x_ref,d,di,...
                 v_feas,mpc);
            end
            x = xhat;

        end
        iter = iter+1;
    end

    % If feas. slack variable is positive, the feasibility solver did not
    % found a feasible point
    if v_feas > 0
        mpc.unfeasible = 1;
    end

end