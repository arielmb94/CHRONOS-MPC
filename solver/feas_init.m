% This function adapts the gradients of the original mpc problem to be used
% on the feasibility optimization problem
%
% Grad(fi_feas) = [Grad_fi_mpc; -1]
%
function feas_slv = feas_init(mpc)

% state inequalities
if ~isempty(mpc.s_cnstr)

    if ~isempty(mpc.s_cnstr.min)
        [b,a] = size(mpc.s_cnstr.grad_min);
        feas_slv.gradXmin = [mpc.s_cnstr.grad_min;-ones(1,a)];
    
        feas_slv.hessXmin = genHessIneq(feas_slv.gradXmin);
    end

    if ~isempty(mpc.s_cnstr.max) 
        [b,a] = size(mpc.s_cnstr.grad_max);
        feas_slv.gradXmax = [mpc.s_cnstr.grad_max;-ones(1,a)];
    
        feas_slv.hessXmax = genHessIneq(feas_slv.gradXmax);
    end
end

% terminal state inequalities
if ~isempty(mpc.s_ter_cnstr)

    if ~isempty(mpc.s_ter_cnstr.min)
        [b,a] = size(mpc.s_ter_cnstr.grad_min);
        feas_slv.gradXtermin = [mpc.s_ter_cnstr.grad_min;-ones(1,a)];
    
        feas_slv.hessXtermin = genHessIneq(feas_slv.gradXtermin);
    end

    if ~isempty(mpc.s_ter_cnstr.max)
        [b,a] = size(mpc.s_ter_cnstr.grad_max);
        feas_slv.gradXtermax = [mpc.s_ter_cnstr.grad_max;-ones(1,a)];
    
        feas_slv.hessXtermax = genHessIneq(feas_slv.gradXtermax);
    end

end

% control inequalities
if ~isempty(mpc.u_cnstr)

    if ~isempty(mpc.u_cnstr.min)
        [b,a] = size(mpc.u_cnstr.grad_min);
        feas_slv.gradUmin = [mpc.u_cnstr.grad_min;-ones(1,a)];
    
        feas_slv.hessUmin = genHessIneq(feas_slv.gradUmin);
    end
    
    if ~isempty(mpc.u_cnstr.max)
        [b,a] = size(mpc.u_cnstr.grad_max);
        feas_slv.gradUmax = [mpc.u_cnstr.grad_max;-ones(1,a)];
    
        feas_slv.hessUmax = genHessIneq(feas_slv.gradUmax);
    end
end

% control differential inequalities
if ~isempty(mpc.du_cnstr)

    if ~isempty(mpc.du_cnstr.min)
        [b,a] = size(mpc.du_cnstr.grad_min);
        feas_slv.gradDeltaUmin = [mpc.du_cnstr.grad_min;-ones(1,a)];
    
        feas_slv.hessDeltaUmin = genHessIneq(feas_slv.gradDeltaUmin);
    end
    
    if ~isempty(mpc.du_cnstr.max)
        [b,a] = size(mpc.du_cnstr.grad_max);
        feas_slv.gradDeltaUmax = [mpc.du_cnstr.grad_max;-ones(1,a)];
    
        feas_slv.hessDeltaUmax = genHessIneq(feas_slv.gradDeltaUmax);
    end
end

% output inequalities
if ~isempty(mpc.y_cnstr)

    if ~isempty(mpc.y_cnstr.min)
        [b,a] = size(mpc.y_cnstr.grad_min);
        feas_slv.gradYmin = [mpc.y_cnstr.grad_min;-ones(1,a)];
    
        feas_slv.hessYmin = genHessIneq(feas_slv.gradYmin);
    end
    
    if ~isempty(mpc.y_cnstr.max)
        [b,a] = size(mpc.y_cnstr.grad_max);
        feas_slv.gradYmax = [mpc.y_cnstr.grad_max;-ones(1,a)];
    
        feas_slv.hessYmax = genHessIneq(feas_slv.gradYmax);
    end
end

% General Linear inequalities
if ~isempty(mpc.h_cnstr)

    if ~isempty(mpc.h_cnstr.min)
        [b,a] = size(mpc.h_cnstr.grad_min);
        feas_slv.gradHmin = [mpc.h_cnstr.grad_min;-ones(1,a)];
    
        feas_slv.hessHmin = genHessIneq(feas_slv.gradHmin);
    end
    
    if ~isempty(mpc.h_cnstr.max)
        [b,a] = size(mpc.h_cnstr.grad_max);
        feas_slv.gradHmax = [mpc.h_cnstr.grad_max;-ones(1,a)];
    
        feas_slv.hessHmax = genHessIneq(feas_slv.gradHmax);
    end
end

end