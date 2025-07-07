%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [u0,x,iter,iter_feas] = mpc_solve(mpc,x0,s_prev,u_prev,r,d,x_ref,dz,dh)
%
% Solve the current iteration of the MPC problem.
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - x0: Nx+Nu column vector, initial guess solution for CHRONOS interior 
%   point iterative solver
%   - s_prev: nx column vector, last measured or estimated system state
%   value
%   - u_prev: nu column vector, control action applied to the system on the
%   previous sampling time
%   - r (optional): tracking reference for the MPC. It can be a ny column 
%   vector (the same reference applies for the full prediction horizon) or 
%   can be a Ny column vector (the user passes a unique reference for each 
%   step of the prediction horizon). If not used, the user must pass an 
%   empty vector [].
%   - d (optional): disturbance input to the system dynamics and to the
%   output signal y models. It can be a nd column vector (the same
%   disturbance applies for the full prediction horizon) or can be an Nd
%   column vector (the user passes a unique disturbance for each step of 
%   the prediction horizon). If not used, the user must pass an empty
%   vector [].
%   - x_ref (optional): nx column vector, reference for the terminal state
%   xN of the prediction horizon. Required whenever the MPC problem
%   contains terminal ingredients. If not used, the user must pass an empty
%   vector [].
%   - dz (optional): disturbance input to the user defined signal model z 
%   for custom cost functions. It can be a ndz column vector (the same
%   disturbance applies for the full prediction horizon) or can be an Ndz
%   column vector (the user passes a unique disturbance for each step of 
%   the prediction horizon). If not used, the user must pass an empty
%   vector [].
%   - dh (optional): disturbance input to the user defined signal model h 
%   for custom constraints. It can be a ndh column vector (the same
%   disturbance applies for the full prediction horizon) or can be an Ndh
%   column vector (the user passes a unique disturbance for each step of 
%   the prediction horizon). If not used, the user must pass an empty
%   vector [].
%
% Out:
%   - u0: nu column vector, first step of the control action sequence
%   computed as solution to the MPC problem
%   - x0: Nx+Nu column vector, optimization variables solution vector to
%   the MPC problem
%   - iter: number of iterations required for the MPC optimization problem 
%   - iter_feas: number of iterations required for the step 0 feasibility
%   starting point finder 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u0,x0,iter,iter_feas] = mpc_solve(mpc,x0,s_prev,u_prev,...
                                            r_in,d_in,x_ref,dz_in,dh_in)

    % number of variables
    n = length(x0);

    % number of equality constraints
    n_eq = size(mpc.Aeq,1); 

    % handle input vector sizes
    if ~isempty(r_in) && length(r_in) < mpc.Ny
        r = fill_vec(r_in,mpc.ny,mpc.Ny,1);
    else
        r = r_in;
    end

    if ~isempty(d_in) && length(d_in)< mpc.Nd
        d = fill_vec(d_in,mpc.nd,mpc.Nd,1);
    else
        d = d_in;
    end

    if ~isempty(dz_in) && length(dz_in)< mpc.Nd
        dz = fill_vec(dz_in,mpc.ndz,mpc.Ndz,1);
    else
        dz = dz_in;
    end

    if ~isempty(dh_in) && length(dh_in)< mpc.Nd
        dh = fill_vec(dh_in,mpc.ndh,mpc.Ndh,1);
    else
        dh = dh_in;
    end

    if mpc.ter_ingredients
        if mpc.x_ref_is_y && isempty(x_ref)
            x_ref = r(end-mpc.ny+1:end);
        end    
    else 
        %x_ref = [];
        grad_ter = [];
    end

    % update b matrix from equality condition
    mpc = update_mpc_beq(mpc,s_prev,d);

    if isempty(mpc.gradPerfQz)
        perfCost = 0;
        z = [];
    else
        perfCost = 1;
    end
    
    
    % Recompute hessian if cost terms have been updated
    if mpc.recompute_cost_hess
        mpc = update_mpc_f0_hess(mpc);
    end

    % Set Newton solver condition at start
    continue_Newton = true;
    iter = 0;
    iter_feas = 0;

    % check if initial state is feasible, it might lead to numerical
    % problems otherwise
    feas = feas_check_s_prev(s_prev,mpc);

    if feas
        % get mpc variables from optimization vector x and constraint
        % information and feasibility
        [mpc,s,s_all,s_ter,u,du,y,err,h,feas] = ...
        get_state_constraint_info(x0,s_prev,u_prev,r,x_ref,d,dh,mpc);
    end

    if ~feas
        % If initial point was not feasible, compute new feasible point
        [x0,iter_feas,mpc] = feas_solve(x0,mpc,s_prev,u_prev,d,x_ref,dh);

        % If feas. solver fails skip mpc solver
        if mpc.unfeasible
            continue_Newton = false;
        else
            iter = iter+iter_feas;
            % Get constraint info on new point
            [mpc,s,s_all,s_ter,u,du,y,err,h,feas] = ...
                get_state_constraint_info(x0,s_prev,u_prev,r,x_ref,d,dh,mpc);
        end

    end

    opts.SYM = true;
    lambda2 = 1;

    % if the number of feasibility solver iterations surpassed the limit
    % for mpc iterations but finds a feasible starting point, reduce the
    % iter counter to at least execute a newton iteration
    if iter > mpc.max_iter
        iter = mpc.max_iter-1;
    end

    while mpc.eps <= lambda2*0.5 && continue_Newton && iter < mpc.max_iter

        % Compute gradient:

        % 1. Compute gradient of box inequalities at x0:
        % init inequalities gradient vector
        grad_fi_Ind = zeros(n,1);

        % state inequalities
        if ~isempty(mpc.s_cnstr)
            if ~isempty(mpc.s_cnstr.min)
                grad_s_min_Ind_x0 = grad_box_Ind(mpc.s_cnstr.fi_min_x0,...
                                                 mpc.s_cnstr.grad_min);
    
                grad_fi_Ind = grad_fi_Ind + grad_s_min_Ind_x0;
            end

            if ~isempty(mpc.s_cnstr.max)
                grad_x_max_Ind_x0 = grad_box_Ind(mpc.s_cnstr.fi_max_x0,...
                                                 mpc.s_cnstr.grad_max);
    
                grad_fi_Ind = grad_fi_Ind + grad_x_max_Ind_x0;
            end
        end

        % terminal state inequalities
        if ~isempty(mpc.s_ter_cnstr)              
            if  ~isempty(mpc.s_ter_cnstr.min)
                grad_s_ter_min_Ind_x0 = grad_box_Ind(mpc.s_ter_cnstr.fi_min_x0,...
                                                    mpc.s_ter_cnstr.grad_min);
    
                grad_fi_Ind = grad_fi_Ind + grad_s_ter_min_Ind_x0;
            end

            if ~isempty(mpc.s_ter_cnstr.max)
                grad_x_ter_max_Ind_x0 = grad_box_Ind(mpc.s_ter_cnstr.fi_max_x0,...
                                                    mpc.s_ter_cnstr.grad_max);
    
                grad_fi_Ind = grad_fi_Ind + grad_x_ter_max_Ind_x0;
            end
        end

        % control inequalities
        if ~isempty(mpc.u_cnstr)
            if ~isempty(mpc.u_cnstr.min)
                grad_u_min_Ind_x0 = grad_box_Ind(mpc.u_cnstr.fi_min_x0,...
                                                 mpc.u_cnstr.grad_min);
    
    
                grad_fi_Ind = grad_fi_Ind + grad_u_min_Ind_x0;
            end
    
            if ~isempty(mpc.u_cnstr.max)
                grad_u_max_Ind_x0 = grad_box_Ind(mpc.u_cnstr.fi_max_x0,...
                                                 mpc.u_cnstr.grad_max);
    
                grad_fi_Ind = grad_fi_Ind + grad_u_max_Ind_x0;
            end
        end

        % control differential inequalities
        if ~isempty(mpc.du_cnstr)
            if ~isempty(mpc.du_cnstr.min)
                grad_du_min_Ind_x0 = grad_box_Ind(mpc.du_cnstr.fi_min_x0,...
                                                 mpc.du_cnstr.grad_min);
    
                grad_fi_Ind = grad_fi_Ind + grad_du_min_Ind_x0;
            end
    
            if ~isempty(mpc.du_cnstr.max)
                grad_du_max_Ind_x0 = grad_box_Ind(mpc.du_cnstr.fi_max_x0,...
                                                 mpc.du_cnstr.grad_max);
    
                grad_fi_Ind = grad_fi_Ind + grad_du_max_Ind_x0;
            end
        end

        % output inequalities
        if ~isempty(mpc.y_cnstr)
            if ~isempty(mpc.y_cnstr.min)
                grad_y_min_Ind_x0 = grad_box_Ind(mpc.y_cnstr.fi_min_x0,...
                                                 mpc.y_cnstr.grad_min);
    
                grad_fi_Ind = grad_fi_Ind + grad_y_min_Ind_x0;
            end
    
            if ~isempty(mpc.y_cnstr.max)
                grad_y_max_Ind_x0 = grad_box_Ind(mpc.y_cnstr.fi_max_x0,...
                                                 mpc.y_cnstr.grad_max);
    
                grad_fi_Ind = grad_fi_Ind + grad_y_max_Ind_x0;
            end
        end

        % General Linear inequalities
        if ~isempty(mpc.h_cnstr)
            if ~isempty(mpc.h_cnstr.min)
                grad_h_min_Ind_x0 = grad_box_Ind(mpc.h_cnstr.fi_min_x0,...
                                                 mpc.h_cnstr.grad_min);
    
                grad_fi_Ind = grad_fi_Ind + grad_h_min_Ind_x0;
            end
    
            if ~isempty(mpc.h_cnstr.max)
                grad_h_max_Ind_x0 = grad_box_Ind(mpc.h_cnstr.fi_max_x0,...
                                                 mpc.h_cnstr.grad_max);
    
                grad_fi_Ind = grad_fi_Ind + grad_h_max_Ind_x0;
            end
        end

        % 2. If enabled, compute terminal ingredients 
        if mpc.ter_ingredients
            [grad_ter,grad_ter_Ind_x0,hess_ter_Ind_x0] = ...
                ter_set_Ind_fun(x_ref,s_ter,mpc.fi_ter_x0,...
                mpc.P,mpc.Nx,mpc.Nu,mpc.nx,mpc.ter_constraint);
            if mpc.ter_constraint
                grad_fi_Ind = grad_fi_Ind + grad_ter_Ind_x0; 
            end
        end
        
        % 3. Compute gradient of cost function at x0
        if perfCost
            z = get_lin_out(s_all,u,dz,mpc.nx,mpc.nu,mpc.nz,mpc.ndz,...
                mpc.N_ctr_hor,mpc.Nz,mpc.Cz,mpc.Dz,mpc.Ddz);
        end

        grad_f0 = grad_f0_MPC(mpc,err,du,u,grad_ter,z);
        
        % 4. Compute gradient at x0 : grad(J) = t*grad(f0)+grad(Phi)
        grad_J_x0 = mpc.t*grad_f0+grad_fi_Ind;


        % Compute Hessian Matrix

        % 1. Compute Hessian of box inequalities at x0:

        % init inequalities hessian vector
        hess_fi_Ind = zeros(n);

        % state inequalities
        if ~isempty(mpc.s_cnstr.min)
            if ~isempty(mpc.s_cnstr.min)
                hess_s_min_Ind_x0 = hess_linear_Ind(mpc.s_cnstr.fi_min_x0,...
                                                    mpc.s_cnstr.hess_min);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_min_Ind_x0;
            end
    
            if ~isempty(mpc.s_cnstr.max)
                hess_s_max_Ind_x0 = hess_linear_Ind(mpc.s_cnstr.fi_max_x0,...
                                                    mpc.s_cnstr.hess_max);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_max_Ind_x0;
            end
        end

        % terminal state inequalities
        if ~isempty(mpc.s_ter_cnstr)
            if ~isempty(mpc.s_ter_cnstr.min)
                hess_s_ter_min_Ind_x0 = hess_linear_Ind(mpc.s_ter_cnstr.fi_min_x0,...
                                                    mpc.s_ter_cnstr.hess_min);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_ter_min_Ind_x0;
            end
    
            if ~isempty(mpc.s_ter_cnstr.max)
                hess_s_ter_max_Ind_x0 = hess_linear_Ind(mpc.s_ter_cnstr.fi_max_x0,...
                                                    mpc.s_ter_cnstr.hess_max);
    
                hess_fi_Ind = hess_fi_Ind + hess_s_ter_max_Ind_x0;
            end
        end

        % control inequalities
        if ~isempty(mpc.u_cnstr)
            if ~isempty(mpc.u_cnstr.min)
                hess_u_min_Ind_x0 = hess_linear_Ind(mpc.u_cnstr.fi_min_x0,...
                                                    mpc.u_cnstr.hess_min);
    
                hess_fi_Ind = hess_fi_Ind + hess_u_min_Ind_x0;
            end
    
            if ~isempty(mpc.u_cnstr.max)
                hess_u_max_Ind_x0 = hess_linear_Ind(mpc.u_cnstr.fi_max_x0,...
                                                    mpc.u_cnstr.hess_max);
    
                hess_fi_Ind = hess_fi_Ind + hess_u_max_Ind_x0;
            end
        end

        % control differential inequalities
        if ~isempty(mpc.du_cnstr)
            if ~isempty(mpc.du_cnstr.min)
                hess_du_min_Ind_x0 = hess_linear_Ind(mpc.du_cnstr.fi_min_x0,...
                                                    mpc.du_cnstr.hess_min);
    
                hess_fi_Ind = hess_fi_Ind + hess_du_min_Ind_x0;
            end
    
            if ~isempty(mpc.du_cnstr.max)
                hess_du_max_Ind_x0 = hess_linear_Ind(mpc.du_cnstr.fi_max_x0,...
                                                    mpc.du_cnstr.hess_max);
    
                hess_fi_Ind = hess_fi_Ind + hess_du_max_Ind_x0;
            end
        end

        % output inequalities
        if ~isempty(mpc.y_cnstr)
            if ~isempty(mpc.y_cnstr.min)
                hess_y_min_Ind_x0 = hess_linear_Ind(mpc.y_cnstr.fi_min_x0,...
                                                    mpc.y_cnstr.hess_min);
    
                hess_fi_Ind = hess_fi_Ind + hess_y_min_Ind_x0;
            end
    
            if ~isempty(mpc.y_cnstr.max)
                hess_y_max_Ind_x0 = hess_linear_Ind(mpc.y_cnstr.fi_max_x0,...
                                                    mpc.y_cnstr.hess_max);
    
                hess_fi_Ind = hess_fi_Ind + hess_y_max_Ind_x0;
            end
        end

        % General Linear inequalities
        if ~isempty(mpc.h_cnstr)
            if ~isempty(mpc.h_cnstr.min)
                hess_h_min_Ind_x0 = hess_linear_Ind(mpc.h_cnstr.fi_min_x0,...
                                                    mpc.h_cnstr.hess_min);
    
                hess_fi_Ind = hess_fi_Ind + hess_h_min_Ind_x0;
            end
    
            if ~isempty(mpc.h_cnstr.max)
                hess_h_max_Ind_x0 = hess_linear_Ind(mpc.h_cnstr.fi_max_x0,...
                                                    mpc.h_cnstr.hess_max);
    
                hess_fi_Ind = hess_fi_Ind + hess_h_max_Ind_x0;
            end
        end

        % 2. If enabled, add terminal constraint hessian term
        if mpc.ter_ingredients && mpc.ter_constraint
            hess_fi_Ind = hess_fi_Ind + hess_ter_Ind_x0;
        end

        % 3. Compute Hessian of f(x0,t):
        hess_J_x0 = mpc.t*mpc.hessCost+hess_fi_Ind;

        % solve KKT system
        KKT = [hess_J_x0 mpc.Aeq';mpc.Aeq zeros(n_eq)];

        delta_x = - linsolve(KKT,[grad_J_x0;mpc.Aeq*x0-mpc.beq],opts);
        %delta_x = - linsolve(KKT,[grad_J_x0;zeros(n_eq,1)],opts);
        %delta_x = - KKT\[grad_J_x0;mpc.Aeq*x-mpc.beq];
        delta_x_prim = delta_x(1:n);

        % compute lambda^2
        lambda2 = -grad_J_x0'*delta_x_prim;

        % Feasibility line search

        l = 1;
        xhat = x0+l*delta_x_prim;

        [mpc,s,s_all,s_ter,u,du,y,err,h,feas] = ...
        get_state_constraint_info(xhat,s_prev,u_prev,r,x_ref,d,dh,mpc);

        if feas
            x0 = xhat;
        else
            while ~feas

                l = l*mpc.Beta;

                xhat = x0+l*delta_x_prim;

                [mpc,s,s_all,s_ter,u,du,y,err,h,feas] = ...
                get_state_constraint_info(xhat,s_prev,u_prev,r,x_ref,d,dh,mpc);
            end
            x0 = xhat;

            if l<mpc.min_l
                continue_Newton = false;
            end
        end
        iter = iter+1;
    end
  
    if mpc.unfeasible
        % if mpc constraints are unfeasible, return first control action 
        % from the optimization vector, reset optimization vector 
        % afterwards
        u0 = x0(1:mpc.nu);
        x0 = x0*0;
    else
        % Get first control action
        u0 = u(1:mpc.nu);
    end

end