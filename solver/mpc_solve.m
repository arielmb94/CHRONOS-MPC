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
%function [u0,x0,iter,iter_feas,mpc] = mpc_solve(mpc,x0,s_prev,u_prev,...
%                                            r_in,d_in,x_ref_in,dz_in,dh_in)
function [u0,x0,iter,iter_feas,mpc] = mpc_solve(mpc,x0,s_prev,u_prev,...
                                                r_in,d_in,x_ref_in,dz_in,dh_in,...
                                                A0,Bd0,Pk,n_rho,varargin)

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
        if mpc.x_ref_is_y && isempty(x_ref_in)
            x_ref = r(end-mpc.ny+1:end);
        else
            x_ref = x_ref_in;
        end    
        grad_ter = zeros(mpc.nx,1);
    else 
        x_ref = [];
        grad_ter = [];
    end

    % update b matrix from equality condition
    if nargin > 9
        mpc = update_mpc_beq(mpc,s_prev,d,A0,Bd0,Pk,n_rho,varargin{:});
    else
        mpc = update_mpc_beq(mpc,s_prev,d);
    end
    
    % Recompute hessian if cost terms have been updated
    if mpc.recompute_cost_hess
        mpc = update_mpc_f0_hess(mpc);
    end

    % Set Newton solver condition at start
    continue_Newton = true;
    iter = 0;

    % get mpc variables from optimization vector x and constraint
    % information and feasibility
    mpc = get_mpc_variables(mpc,x0,s_prev,u_prev,r,d,dh,dz);
    % If exists, update slack variables
    if mpc.Nv
        [mpc,x0] = mpc_slack_update(mpc,x0,x_ref);
    end
    % compute fi
    [mpc,~] = check_mpc_feasibility(mpc,x_ref);
    
    opts.SYM = true;
    lambda2 = 1;

    while mpc.eps <= lambda2*0.5 && continue_Newton && iter < mpc.max_iter

        % Compute gradient:

        % 1. Compute gradient/Hessian of box inequalities at x0:
        % init inequalities gradient vector
        grad_fi_Ind = zeros(n,1);
        % init inequalities hessian vector
        hess_fi_Ind = zeros(n);

        % state inequalities
        if ~isempty(mpc.s_cnstr)
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_box_Indicator_fun(mpc.s_cnstr,grad_fi_Ind,hess_fi_Ind);
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_slack_Indicator_fun(mpc.s_cnstr,grad_fi_Ind,hess_fi_Ind);
        end

        % control inequalities
        if ~isempty(mpc.u_cnstr)
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_box_Indicator_fun(mpc.u_cnstr,grad_fi_Ind,hess_fi_Ind);
        end

        % control differential inequalities
        if ~isempty(mpc.du_cnstr)
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_box_Indicator_fun(mpc.du_cnstr,grad_fi_Ind,hess_fi_Ind);
        end

        % output inequalities
        if ~isempty(mpc.y_cnstr)
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_box_Indicator_fun(mpc.y_cnstr,grad_fi_Ind,hess_fi_Ind);
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_slack_Indicator_fun(mpc.y_cnstr,grad_fi_Ind,hess_fi_Ind);
        end

        % General Linear inequalities
        if ~isempty(mpc.h_cnstr)
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_box_Indicator_fun(mpc.h_cnstr,grad_fi_Ind,hess_fi_Ind);
            [grad_fi_Ind(:),hess_fi_Ind(:,:)] = ...
                gradient_Hessian_slack_Indicator_fun(mpc.h_cnstr,grad_fi_Ind,hess_fi_Ind);
        end

        % 2. If enabled, compute terminal ingredients 
        % CODEGEN NOTE: to be commented out if there is not terminal ingredients
        if mpc.ter_ingredients
            grad_ter = ter_set_Ind_fun(mpc,x_ref);
%             if mpc.ter_constraint
%                 % add gradient from terminal constraint
%                 grad_fi_Ind = grad_fi_Ind + grad_ter_Ind_x0;
% 
%                 % add gradient from terminal constraint slack variable
%                 % positivity constraint
%                 grad_ter_slack_Ind_x0 = grad_box_Ind(mpc.fi_ter_slack_positivity_x0,...
%                                             mpc.fi_ter_slack_positivity_grad);
% 
%                 grad_fi_Ind = grad_fi_Ind + grad_ter_slack_Ind_x0;
% 
%                 hess_fi_Ind = hess_fi_Ind + hess_ter_Ind_x0;
% 
%                 hess_ter_slack_Ind_x0 = hess_linear_Ind(mpc.fi_ter_slack_positivity_x0,...
%                                             mpc.fi_ter_slack_positivity_hess);
% 
%                 hess_fi_Ind = hess_fi_Ind + hess_ter_slack_Ind_x0;
%             end
        end

        grad_f0 = grad_f0_MPC(mpc,mpc.err,mpc.du,mpc.u,grad_ter,mpc.z);
        
        % 4. Compute gradient at x0 : grad(J) = t*grad(f0)+grad(Phi)
        grad_J_x0 = mpc.t*grad_f0+grad_fi_Ind;


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

        mpc = get_mpc_variables(mpc,xhat,s_prev,u_prev,r,d,dh,dz);
        [mpc,feas] = check_mpc_feasibility(mpc,x_ref);

        if feas
            x0 = xhat;
        else
            while ~feas
                l = l*mpc.Beta;

                xhat = x0+l*delta_x_prim;

                mpc = get_mpc_variables(mpc,xhat,s_prev,u_prev,r,d,dh,dz);
                [mpc,feas] = check_mpc_feasibility(mpc,x_ref);

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
        u0 = mpc.u(1:mpc.nu);
    end

end