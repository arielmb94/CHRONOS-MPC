% This function adapts the gradients of the original mpc problem to be used
% on the feasibility optimization problem
%
% Grad(fi_feas) = [Grad(fi_mpc; -1]
%
function feas_slv = feas_init(mpc)

% state inequalities
if ~isempty(mpc.x_min)
    [b,a] = size(mpc.gradXmin);
    feas_slv.gradXmin = [mpc.gradXmin;-ones(1,a)];

    feas_slv.hessXmin = genHessIneq(feas_slv.gradXmin);
end

if ~isempty(mpc.x_max) 
    [b,a] = size(mpc.gradXmax);
    feas_slv.gradXmax = [mpc.gradXmax;-ones(1,a)];

    feas_slv.hessXmax = genHessIneq(feas_slv.gradXmax);
end

% terminal state inequalities
if ~isempty(mpc.x_ter_min)
    [b,a] = size(mpc.gradXtermin);
    feas_slv.gradXtermin = [mpc.gradXtermin;-ones(1,a)];

    feas_slv.hessXtermin = genHessIneq(feas_slv.gradXtermin);
end

if ~isempty(mpc.x_ter_max)
    [b,a] = size(mpc.gradXtermax);
    feas_slv.gradXtermax = [mpc.gradXtermax;-ones(1,a)];

    feas_slv.hessXtermax = genHessIneq(feas_slv.gradXtermax);
end

% control inequalities
if ~isempty(mpc.u_min)
    [b,a] = size(mpc.gradUmin);
    feas_slv.gradUmin = [mpc.gradUmin;-ones(1,a)];

    feas_slv.hessUmin = genHessIneq(feas_slv.gradUmin);
end

if ~isempty(mpc.u_max)
    [b,a] = size(mpc.gradUmax);
    feas_slv.gradUmax = [mpc.gradUmax;-ones(1,a)];

    feas_slv.hessUmax = genHessIneq(feas_slv.gradUmax);
end

% control differential inequalities
if ~isempty(mpc.du_min)
    [b,a] = size(mpc.gradDeltaUmin);
    feas_slv.gradDeltaUmin = [mpc.gradDeltaUmin;-ones(1,a)];

    feas_slv.hessDeltaUmin = genHessIneq(feas_slv.gradDeltaUmin);
end

if ~isempty(mpc.du_max)
    [b,a] = size(mpc.gradDeltaUmax);
    feas_slv.gradDeltaUmax = [mpc.gradDeltaUmax;-ones(1,a)];

    feas_slv.hessDeltaUmax = genHessIneq(feas_slv.gradDeltaUmax);
end

% output inequalities
if ~isempty(mpc.y_min)
    [b,a] = size(mpc.gradYmin);
    feas_slv.gradYmin = [mpc.gradYmin;-ones(1,a)];

    feas_slv.hessYmin = genHessIneq(feas_slv.gradYmin);
end

if ~isempty(mpc.y_max)
    [b,a] = size(mpc.gradYmax);
    feas_slv.gradYmax = [mpc.gradYmax;-ones(1,a)];

    feas_slv.hessYmax = genHessIneq(feas_slv.gradYmax);
end

% General Linear inequalities
if ~isempty(mpc.h_min)
    [b,a] = size(mpc.gradHmin);
    feas_slv.gradHmin = [mpc.gradHmin;-ones(1,a)];

    feas_slv.hessHmin = genHessIneq(feas_slv.gradHmin);
end

if ~isempty(mpc.h_max)
    [b,a] = size(mpc.gradHmax);
    feas_slv.gradHmax = [mpc.gradHmax;-ones(1,a)];

    feas_slv.hessHmax = genHessIneq(feas_slv.gradHmax);
end

end