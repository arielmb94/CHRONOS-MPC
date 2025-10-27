% This function adapts the gradients of the original mpc problem to be used
% on the feasibility optimization problem
%
% Grad(fi_feas) = [Grad_fi_mpc; -1]
%
function mpc = feas_init(mpc)

% output inequalities
if ~isempty(mpc.y_cnstr)

    if ~isempty(mpc.y_cnstr.min)
        a = size(mpc.y_cnstr.grad_min,2);
        mpc.y_cnstr.grad_min_feas_slv(:)  = [mpc.y_cnstr.grad_min;-ones(1,a)];
    
        mpc.y_cnstr.hess_min_feas_slv = updateHessIneq(mpc.y_cnstr.grad_min_feas_slv, ...
                                            mpc.y_cnstr.hess_min_feas_slv);
    end
    
    if ~isempty(mpc.y_cnstr.max)
        a = size(mpc.y_cnstr.grad_max,2);
        mpc.y_cnstr.grad_max_feas_slv(:)  = [mpc.y_cnstr.grad_max;-ones(1,a)];
    
        mpc.y_cnstr.hess_max_feas_slv = updateHessIneq(mpc.y_cnstr.grad_max_feas_slv, ...
                                            mpc.y_cnstr.hess_max_feas_slv);
    end
end

% General Linear inequalities
if ~isempty(mpc.h_cnstr)

    if ~isempty(mpc.h_cnstr.min)
        a = size(mpc.h_cnstr.grad_min,2);
        mpc.h_cnstr.grad_min_feas_slv(:)  = [mpc.h_cnstr.grad_min;-ones(1,a)];
    
        mpc.h_cnstr.hess_min_feas_slv = updateHessIneq(mpc.h_cnstr.grad_min_feas_slv, ...
                                            mpc.h_cnstr.hess_min_feas_slv);
    end
    
    if ~isempty(mpc.h_cnstr.max)
        a = size(mpc.h_cnstr.grad_max,2);
        mpc.h_cnstr.grad_max_feas_slv(:)  = [mpc.h_cnstr.grad_max;-ones(1,a)];
    
        mpc.h_cnstr.hess_max_feas_slv = updateHessIneq(mpc.h_cnstr.grad_max_feas_slv, ...
                                            mpc.h_cnstr.hess_max_feas_slv);
    end
end

end