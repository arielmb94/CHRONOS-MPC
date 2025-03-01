function [u0,x] = feas_solve(x0,s_prev,u_prev,r,d,mpc,x_ref,di)

    % extend vector of optimization variables with feasibility slack
    x = [x0;mpc.x0_feas];

    % number of variables
    n = length(x);
    
    % feasibility cost function is J = v
    % define fesibility cost function gradient as [0 ... 0 1]^T
    grad_f0 = [zeros(n-1,1);1];

    % number of equality constraints
    n_eq = size(mpc.Aeq,1); 

    % update b matrix from equality condition
    mpc = update_mpc_beq(mpc,s_prev,d);

    if mpc.ter_ingredients
        if mpc.x_ref_is_y && isempty(x_ref)
            x_ref = r;
            % TODO: account for reference being a sequence
        end    
    else 
        x_ref = [];
        grad_ter = [];
    end

    % for first iteration we assume x0 is feasible and wont check
    check_feas = false;
    % get mpc variables from optimimization vector x and constraint
    % information
    [s,s_all,s_ter,u,du,y,yi,...
    fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
    fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
    fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
    fi_yi_min_x0,fi_yi_max_x0,feas] = ...
    get_feas_probl_constraint_info(x,s_prev,u_prev,x_ref,d,di,mpc,check_feas);
    % for following iteration check feasibility
    check_feas = true;

    continue_Newton = true;
    opts.SYM = true;
    lambda2 = 1;
    while mpc.eps <= lambda2*0.5 && continue_Newton

        % Compute gradient:

        % 1. Compute gradient of box inequalities at x0:
        % init inequalities gradient vector
        grad_fi_Ind = zeros(n,1);

        % state inequalities
        if ~isempty(mpc.x_min)
            grad_s_min_Ind_x0 = grad_box_Ind(fi_s_min_x0,mpc.gradXmin);

            grad_fi_Ind = grad_fi_Ind + grad_s_min_Ind_x0;
        end

        if ~isempty(mpc.x_max)
            grad_x_max_Ind_x0 = grad_box_Ind(fi_s_max_x0,mpc.gradXmax);

            grad_fi_Ind = grad_fi_Ind + grad_x_max_Ind_x0;
        end

        % terminal state inequalities
        if ~isempty(mpc.x_ter_min)
            grad_s_ter_min_Ind_x0 = grad_box_Ind(fi_s_ter_min_x0,mpc.gradXtermin);

            grad_fi_Ind = grad_fi_Ind + grad_s_ter_min_Ind_x0;
        end

        if ~isempty(mpc.x_ter_max)
            grad_x_ter_max_Ind_x0 = grad_box_Ind(fi_s_ter_max_x0,mpc.gradXtermax);

            grad_fi_Ind = grad_fi_Ind + grad_x_ter_max_Ind_x0;
        end

        % control inequalities
        if ~isempty(mpc.u_min)
            grad_u_min_Ind_x0 = grad_box_Ind(fi_u_min_x0,mpc.gradUmin);

            grad_fi_Ind = grad_fi_Ind + grad_u_min_Ind_x0;
        end

        if ~isempty(mpc.u_max)
            grad_u_max_Ind_x0 = grad_box_Ind(fi_u_max_x0,mpc.gradUmax);

            grad_fi_Ind = grad_fi_Ind + grad_u_max_Ind_x0;
        end

        % control differential inequalities
        if ~isempty(mpc.du_min)
            grad_du_min_Ind_x0 = grad_box_Ind(fi_du_min_x0,mpc.gradDeltaUmin);

            grad_fi_Ind = grad_fi_Ind + grad_du_min_Ind_x0;
        end

        if ~isempty(mpc.du_max)
            grad_du_max_Ind_x0 = grad_box_Ind(fi_du_max_x0,mpc.gradDeltaUmax);

            grad_fi_Ind = grad_fi_Ind + grad_du_max_Ind_x0;
        end

        % output inequalities
        if ~isempty(mpc.y_min)
            grad_y_min_Ind_x0 = grad_box_Ind(fi_y_min_x0,mpc.gradYmin);

            grad_fi_Ind = grad_fi_Ind + grad_y_min_Ind_x0;
        end

        if ~isempty(mpc.y_max)
            grad_y_max_Ind_x0 = grad_box_Ind(fi_y_max_x0,mpc.gradYmax);

            grad_fi_Ind = grad_fi_Ind + grad_y_max_Ind_x0;
        end

        % General Linear inequalities
        if ~isempty(mpc.yi_min)
            grad_yi_min_Ind_x0 = grad_box_Ind(fi_yi_min_x0,mpc.gradYimin);

            grad_fi_Ind = grad_fi_Ind + grad_yi_min_Ind_x0;
        end

        if ~isempty(mpc.yi_max)
            grad_yi_max_Ind_x0 = grad_box_Ind(fi_yi_max_x0,mpc.gradYimax);

            grad_fi_Ind = grad_fi_Ind + grad_yi_max_Ind_x0;
        end

        % 2. If enabled, compute terminal ingredients 
        if mpc.ter_ingredients
            [grad_ter,grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
                ter_set_Ind_fun(x_ref,s_ter,fi_ter_x0,...
                mpc.P,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.N,mpc.ter_constraint);
            if mpc.ter_constraint
                grad_fi_Ind = grad_fi_Ind + grad_ter_Ind_x0; 
            end
        end
        
        % 4. Compute gradient at x0 : grad(J) = t*grad(f0)+grad(Phi)
        grad_J_x0 = mpc.t_feas*grad_f0+grad_fi_Ind;


        % Compute Hessian Matrix

        % 1. Compute Hessian of box inequalities at x0:

        % init inequalities hessian vector
        hess_fi_Ind = zeros(n);

        % state inequalities
        if ~isempty(mpc.x_min)
            hess_s_min_Ind_x0 = hess_linear_Ind(fi_s_min_x0,mpc.hessXmin);

            hess_fi_Ind = hess_fi_Ind + hess_s_min_Ind_x0;
        end

        if ~isempty(mpc.x_max)
            hess_s_max_Ind_x0 = hess_linear_Ind(fi_s_max_x0,mpc.hessXmax);

            hess_fi_Ind = hess_fi_Ind + hess_s_max_Ind_x0;
        end

        % terminal state inequalities
        if ~isempty(mpc.x_ter_min)
            hess_s_ter_min_Ind_x0 = hess_linear_Ind(fi_s_ter_min_x0,mpc.hessXtermin);

            hess_fi_Ind = hess_fi_Ind + hess_s_ter_min_Ind_x0;
        end

        if ~isempty(mpc.x_ter_max)
            hess_s_ter_max_Ind_x0 = hess_linear_Ind(fi_s_ter_max_x0,mpc.hessXtermax);

            hess_fi_Ind = hess_fi_Ind + hess_s_ter_max_Ind_x0;
        end

        % control inequalities
        if ~isempty(mpc.u_min)
            hess_u_min_Ind_x0 = hess_linear_Ind(fi_u_min_x0,mpc.hessUmin);

            hess_fi_Ind = hess_fi_Ind + hess_u_min_Ind_x0;
        end

        if ~isempty(mpc.u_max)
            hess_u_max_Ind_x0 = hess_linear_Ind(fi_u_max_x0,mpc.hessUmax);

            hess_fi_Ind = hess_fi_Ind + hess_u_max_Ind_x0;
        end

        % control differential inequalities
        if ~isempty(mpc.du_min)
            hess_du_min_Ind_x0 = hess_linear_Ind(fi_du_min_x0,mpc.hessDeltaUmin);

            hess_fi_Ind = hess_fi_Ind + hess_du_min_Ind_x0;
        end

        if ~isempty(mpc.du_max)
            hess_du_max_Ind_x0 = hess_linear_Ind(fi_du_max_x0,mpc.hessDeltaUmax);

            hess_fi_Ind = hess_fi_Ind + hess_du_max_Ind_x0;
        end

        % output inequalities
        if ~isempty(mpc.y_min)
            hess_y_min_Ind_x0 = hess_linear_Ind(fi_y_min_x0,mpc.hessYmin);

            hess_fi_Ind = hess_fi_Ind + hess_y_min_Ind_x0;
        end

        if ~isempty(mpc.y_max)
            hess_y_max_Ind_x0 = hess_linear_Ind(fi_y_max_x0,mpc.hessYmax);

            hess_fi_Ind = hess_fi_Ind + hess_y_max_Ind_x0;
        end

        % General Linear inequalities
        if ~isempty(mpc.yi_min)
            hess_yi_min_Ind_x0 = hess_linear_Ind(fi_yi_min_x0,mpc.hessYimin);

            hess_fi_Ind = hess_fi_Ind + hess_yi_min_Ind_x0;
        end

        if ~isempty(mpc.yi_max)
            hess_yi_max_Ind_x0 = hess_linear_Ind(fi_yi_max_x0,mpc.hessYimax);

            hess_fi_Ind = hess_fi_Ind + hess_yi_max_Ind_x0;
        end

        % 2. If enabled, add terminal constraint hessian term
        if mpc.ter_ingredients && mpc.ter_constraint
            hess_fi_Ind = hess_fi_Ind + hess_ter_Ind_x0;
        end

        % solve KKT system
        KKT = [hess_fi_Ind mpc.Aeq';mpc.Aeq zeros(n_eq)];

        delta_x = - linsolve(KKT,[grad_J_x0;mpc.Aeq*x-mpc.beq],opts);
        %delta_x = - linsolve(KKT,[grad_J_x0;zeros(n_eq,1)],opts);
        %delta_x = - KKT\[grad_J_x0;mpc.Aeq*x-mpc.beq];
        delta_x_prim = delta_x(1:n);

        % compute lambda^2
        lambda2 = -grad_J_x0'*delta_x_prim;

        % Feasibility line search

        l = 1;
        xhat = x+l*delta_x_prim;

        [s,s_all,s_ter,u,du,y,yi,...
         fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
         fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
         fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
         fi_yi_min_x0,fi_yi_max_x0,feas] = ...
         get_feas_probl_constraint_info(xhat,s_prev,u_prev,x_ref,d,di,mpc,check_feas);

        if feas
            x = xhat;
        else
            while ~feas

                l = l*mpc.Beta;

                xhat = x+l*delta_x_prim;

                [s,s_all,s_ter,u,du,y,yi,...
                 fi_s_min_x0,fi_s_max_x0,fi_s_ter_min_x0,fi_s_ter_max_x0,...
                 fi_u_min_x0,fi_u_max_x0,fi_du_min_x0,fi_du_max_x0,...
                 fi_y_min_x0,fi_y_max_x0,fi_ter_x0,...
                 fi_yi_min_x0,fi_yi_max_x0,feas] = ...
                 get_feas_probl_constraint_info(xhat,s_prev,u_prev,x_ref,d,di,mpc,check_feas);
            end
            x = xhat;

            if l<mpc.min_l
                continue_Newton = false;
            end
        end

    end

    % Get first control action
    u0 = u(1:mpc.nu);

    %J = f0_fun_MPC(mpc.Qe,err,mpc.N,mpc.ny,mpc.Rdu,du,mpc.nu,[],[]);
    J = 0;

end