function [x_mpc,iter,mpc] = feas_solve(x0,mpc,s_prev,u_prev,d,x_ref,dh)

    % Set false at start
    mpc.unfeasible = 0;

    x_mpc = x0;
    v_feas = mpc.v0_feas;
    % extend optimization vector with feasibility slack variable
    x = [x_mpc;v_feas];

    % number of variables of original mpc problem
    n = length(x);
    
    % Generate constraint's gradient from original mpc
    mpc = feas_init(mpc);

    % Adapt equality constraints to feasibility solver size 
    % number of equality constraints
    n_eq = size(mpc.Aeq,1); 
    Aeq_feas = [mpc.Aeq zeros(n_eq,1)];

    % get mpc variables from optimization vector x and constraint
    % information
    % increase v_feas if needed   
    feas = 0;
    while ~feas

        [mpc,feas] = get_feas_probl_constraint_info(...
                                            x_mpc,s_prev,u_prev,x_ref,d,dh,...
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
        if ~isempty(mpc.s_cnstr)
            if ~isempty(mpc.s_cnstr.min)
                grad_s_min_Ind_x0 = grad_box_Ind(mpc.s_cnstr.fi_min_x0,...
                                                 mpc.s_cnstr.grad_min_feas_slv,...
                                                 mpc.s_cnstr.min_activ_set,...
                                                 mpc.s_cnstr.min_activ_indicator, ...
                                                 mpc.nx);
    
                grad_fi_Ind = grad_fi_Ind + grad_s_min_Ind_x0;
            end

            if ~isempty(mpc.s_cnstr.max)
                grad_x_max_Ind_x0 = grad_box_Ind(mpc.s_cnstr.fi_max_x0,...
                                                 mpc.s_cnstr.grad_max_feas_slv,...
                                                 mpc.s_cnstr.max_activ_set,...
                                                 mpc.s_cnstr.max_activ_indicator, ...
                                                 mpc.nx);
    
                grad_fi_Ind = grad_fi_Ind + grad_x_max_Ind_x0;
            end
        end

        % terminal state inequalities
        if ~isempty(mpc.s_ter_cnstr)              
            if  ~isempty(mpc.s_ter_cnstr.min)
                grad_s_ter_min_Ind_x0 = grad_box_Ind(mpc.s_ter_cnstr.fi_min_x0,...
                                                    mpc.s_ter_cnstr.grad_min_feas_slv,...
                                                    mpc.s_ter_cnstr.min_activ_set,...
                                                    mpc.s_ter_cnstr.min_activ_indicator, ...
                                                    mpc.nx);
    
                grad_fi_Ind = grad_fi_Ind + grad_s_ter_min_Ind_x0;
            end

            if ~isempty(mpc.s_ter_cnstr.max)
                grad_x_ter_max_Ind_x0 = grad_box_Ind(mpc.s_ter_cnstr.fi_max_x0,...
                                                    mpc.s_ter_cnstr.grad_max_feas_slv,...
                                                    mpc.s_ter_cnstr.max_activ_set,...
                                                    mpc.s_ter_cnstr.max_activ_indicator, ...
                                                    mpc.nx);
    
                grad_fi_Ind = grad_fi_Ind + grad_x_ter_max_Ind_x0;
            end
        end

        % control inequalities
        if ~isempty(mpc.u_cnstr)
            if ~isempty(mpc.u_cnstr.min)
                grad_u_min_Ind_x0 = grad_box_Ind(mpc.u_cnstr.fi_min_x0,...
                                                 mpc.u_cnstr.grad_min_feas_slv,...
                                                 mpc.u_cnstr.min_activ_set,...
                                                 mpc.u_cnstr.min_activ_indicator, ...
                                                 mpc.nu);
    
                grad_fi_Ind = grad_fi_Ind + grad_u_min_Ind_x0;
            end
    
            if ~isempty(mpc.u_cnstr.max)
                grad_u_max_Ind_x0 = grad_box_Ind(mpc.u_cnstr.fi_max_x0,...
                                                 mpc.u_cnstr.grad_max_feas_slv,...
                                                 mpc.u_cnstr.max_activ_set,...
                                                 mpc.u_cnstr.max_activ_indicator, ...
                                                 mpc.nu);
    
                grad_fi_Ind = grad_fi_Ind + grad_u_max_Ind_x0;
            end
        end

        % control differential inequalities
        if ~isempty(mpc.du_cnstr)
            if ~isempty(mpc.du_cnstr.min)
                grad_du_min_Ind_x0 = grad_box_Ind(mpc.du_cnstr.fi_min_x0,...
                                                 mpc.du_cnstr.grad_min_feas_slv,...
                                                 mpc.du_cnstr.min_activ_set,...
                                                 mpc.du_cnstr.min_activ_indicator, ...
                                                 mpc.nu);
    
                grad_fi_Ind = grad_fi_Ind + grad_du_min_Ind_x0;
            end
    
            if ~isempty(mpc.du_cnstr.max)
                grad_du_max_Ind_x0 = grad_box_Ind(mpc.du_cnstr.fi_max_x0,...
                                                 mpc.du_cnstr.grad_max_feas_slv,...
                                                 mpc.du_cnstr.max_activ_set,...
                                                 mpc.du_cnstr.max_activ_indicator, ...
                                                 mpc.nu);
    
                grad_fi_Ind = grad_fi_Ind + grad_du_max_Ind_x0;
            end
        end

        % output inequalities
        if ~isempty(mpc.y_cnstr)
            if ~isempty(mpc.y_cnstr.min)
                grad_y_min_Ind_x0 = grad_box_Ind(mpc.y_cnstr.fi_min_x0,...
                                                 mpc.y_cnstr.grad_min_feas_slv,...
                                                 mpc.y_cnstr.min_activ_set,...
                                                 mpc.y_cnstr.min_activ_indicator, ...
                                                 mpc.ny);
    
                grad_fi_Ind = grad_fi_Ind + grad_y_min_Ind_x0;
            end
    
            if ~isempty(mpc.y_cnstr.max)
                grad_y_max_Ind_x0 = grad_box_Ind(mpc.y_cnstr.fi_max_x0,...
                                                 mpc.y_cnstr.grad_max_feas_slv,...
                                                 mpc.y_cnstr.max_activ_set,...
                                                 mpc.y_cnstr.max_activ_indicator, ...
                                                 mpc.ny);
    
                grad_fi_Ind = grad_fi_Ind + grad_y_max_Ind_x0;
            end
        end

        % General Linear inequalities
        if ~isempty(mpc.h_cnstr)
            if ~isempty(mpc.h_cnstr.min)
                grad_h_min_Ind_x0 = grad_box_Ind(mpc.h_cnstr.fi_min_x0,...
                                                 mpc.h_cnstr.grad_min_feas_slv,...
                                                 mpc.h_cnstr.min_activ_set,...
                                                 mpc.h_cnstr.min_activ_indicator, ...
                                                 mpc.nh);
    
                grad_fi_Ind = grad_fi_Ind + grad_h_min_Ind_x0;
            end
    
            if ~isempty(mpc.h_cnstr.max)
                grad_h_max_Ind_x0 = grad_box_Ind(mpc.h_cnstr.fi_max_x0,...
                                                 mpc.h_cnstr.grad_max_feas_slv,...
                                                 mpc.h_cnstr.max_activ_set,...
                                                 mpc.h_cnstr.max_activ_indicator, ...
                                                 mpc.nh);
    
                grad_fi_Ind = grad_fi_Ind + grad_h_max_Ind_x0;
            end
        end

        % 2. If enabled, compute terminal constraint gradients 
        % CODEGEN NOTE: to be commented out if there is not terminal ingredients
        if mpc.ter_ingredients && mpc.ter_constraint

            [grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
                ter_set_feas_Ind_fun(x_ref,mpc.s_ter,mpc.fi_ter_x0,mpc.P,...
                mpc.nx,n);

                grad_fi_Ind = grad_fi_Ind + grad_ter_Ind_x0; 
        else
            hess_ter_Ind_x0 = zeros(n);
        end


        % Compute Hessian Matrix

        % 1. Compute Hessian of box inequalities at x0:

        % init inequalities hessian vector
        hess_fi_Ind = zeros(n);

        % state inequalities
        if ~isempty(mpc.s_cnstr)
            if ~isempty(mpc.s_cnstr.min)
                hess_s_min_Ind_x0 = hess_linear_Ind(mpc.s_cnstr.fi_min_x0,...
                                                    mpc.s_cnstr.hess_min_feas_slv,...
                                                    mpc.s_cnstr.min_activ_set,...
                                                    mpc.s_cnstr.min_activ_indicator, ...
                                                    mpc.nx);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_min_Ind_x0;
            end
    
            if ~isempty(mpc.s_cnstr.max)
                hess_s_max_Ind_x0 = hess_linear_Ind(mpc.s_cnstr.fi_max_x0,...
                                                    mpc.s_cnstr.hess_max_feas_slv,...
                                                    mpc.s_cnstr.max_activ_set,...
                                                    mpc.s_cnstr.max_activ_indicator, ...
                                                    mpc.nx);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_max_Ind_x0;
            end
        end

        % terminal state inequalities
        if ~isempty(mpc.s_ter_cnstr)
            if ~isempty(mpc.s_ter_cnstr.min)
                hess_s_ter_min_Ind_x0 = hess_linear_Ind(mpc.s_ter_cnstr.fi_min_x0,...
                                                    mpc.s_ter_cnstr.hess_min_feas_slv,...
                                                    mpc.s_ter_cnstr.min_activ_set,...
                                                    mpc.s_ter_cnstr.min_activ_indicator, ...
                                                    mpc.nx);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_ter_min_Ind_x0;
            end
    
            if ~isempty(mpc.s_ter_cnstr.max)
                hess_s_ter_max_Ind_x0 = hess_linear_Ind(mpc.s_ter_cnstr.fi_max_x0,...
                                                    mpc.s_ter_cnstr.hess_max_feas_slv,...
                                                    mpc.s_ter_cnstr.max_activ_set,...
                                                    mpc.s_ter_cnstr.max_activ_indicator, ...
                                                    mpc.nx);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_ter_max_Ind_x0;
            end
        end

        % control inequalities
        if ~isempty(mpc.u_cnstr)
            if ~isempty(mpc.u_cnstr.min)
                hess_u_min_Ind_x0 = hess_linear_Ind(mpc.u_cnstr.fi_min_x0,...
                                                    mpc.u_cnstr.hess_min_feas_slv,...
                                                    mpc.u_cnstr.min_activ_set,...
                                                    mpc.u_cnstr.min_activ_indicator, ...
                                                    mpc.nu);
    
                hess_fi_Ind = hess_fi_Ind + hess_u_min_Ind_x0;
            end
    
            if ~isempty(mpc.u_cnstr.max)
                hess_u_max_Ind_x0 = hess_linear_Ind(mpc.u_cnstr.fi_max_x0,...
                                                    mpc.u_cnstr.hess_max_feas_slv,...
                                                    mpc.u_cnstr.max_activ_set,...
                                                    mpc.u_cnstr.max_activ_indicator, ...
                                                    mpc.nu);
    
                hess_fi_Ind = hess_fi_Ind + hess_u_max_Ind_x0;
            end
        end

        % control differential inequalities
        if ~isempty(mpc.du_cnstr)
            if ~isempty(mpc.du_cnstr.min)
                hess_du_min_Ind_x0 = hess_linear_Ind(mpc.du_cnstr.fi_min_x0,...
                                                    mpc.du_cnstr.hess_min_feas_slv,...
                                                    mpc.du_cnstr.min_activ_set,...
                                                    mpc.du_cnstr.min_activ_indicator, ...
                                                    mpc.nu);
    
                hess_fi_Ind = hess_fi_Ind + hess_du_min_Ind_x0;
            end
    
            if ~isempty(mpc.du_cnstr.max)
                hess_du_max_Ind_x0 = hess_linear_Ind(mpc.du_cnstr.fi_max_x0,...
                                                    mpc.du_cnstr.hess_max_feas_slv,...
                                                    mpc.du_cnstr.max_activ_set,...
                                                    mpc.du_cnstr.max_activ_indicator, ...
                                                    mpc.nu);
    
                hess_fi_Ind = hess_fi_Ind + hess_du_max_Ind_x0;
            end
        end

        % output inequalities
        if ~isempty(mpc.y_cnstr)
            if ~isempty(mpc.y_cnstr.min)
                hess_y_min_Ind_x0 = hess_linear_Ind(mpc.y_cnstr.fi_min_x0,...
                                                    mpc.y_cnstr.hess_min_feas_slv,...
                                                    mpc.y_cnstr.min_activ_set,...
                                                    mpc.y_cnstr.min_activ_indicator, ...
                                                    mpc.ny);
    
                hess_fi_Ind = hess_fi_Ind + hess_y_min_Ind_x0;
            end
    
            if ~isempty(mpc.y_cnstr.max)
                hess_y_max_Ind_x0 = hess_linear_Ind(mpc.y_cnstr.fi_max_x0,...
                                                    mpc.y_cnstr.hess_max_feas_slv,...
                                                    mpc.y_cnstr.max_activ_set,...
                                                    mpc.y_cnstr.max_activ_indicator, ...
                                                    mpc.ny);
    
                hess_fi_Ind = hess_fi_Ind + hess_y_max_Ind_x0;
            end
        end

        % General Linear inequalities
        if ~isempty(mpc.h_cnstr)
            if ~isempty(mpc.h_cnstr.min)
                hess_h_min_Ind_x0 = hess_linear_Ind(mpc.h_cnstr.fi_min_x0,...
                                                    mpc.h_cnstr.hess_min_feas_slv,...
                                                    mpc.h_cnstr.min_activ_set,...
                                                    mpc.h_cnstr.min_activ_indicator, ...
                                                    mpc.nh);
    
                hess_fi_Ind = hess_fi_Ind + hess_h_min_Ind_x0;
            end
    
            if ~isempty(mpc.h_cnstr.max)
                hess_h_max_Ind_x0 = hess_linear_Ind(mpc.h_cnstr.fi_max_x0,...
                                                    mpc.h_cnstr.hess_max_feas_slv,...
                                                    mpc.h_cnstr.max_activ_set,...
                                                    mpc.h_cnstr.max_activ_indicator, ...
                                                    mpc.nh);
    
                hess_fi_Ind = hess_fi_Ind + hess_h_max_Ind_x0;
            end
        end

        % 2. If enabled, add terminal constraint hessian term
        % CODEGEN NOTE: to be commented out if there is not terminal ingredients
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

        [mpc,feas] = get_feas_probl_constraint_info(...
                                            x_mpc,s_prev,u_prev,x_ref,d,dh,...
                                            v_feas,mpc);

        if feas
            x = xhat;
        else
            while ~feas
                
                l = l*mpc.Beta;

                xhat = x+l*delta_x_prim;
                x_mpc = xhat(1:n-1);
                v_feas = xhat(n);

                [mpc,feas] = get_feas_probl_constraint_info(...
                                            x_mpc,s_prev,u_prev,x_ref,d,dh,...
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