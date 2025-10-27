%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_output_cnstr(mpc,y_min,y_max)
%
% Define constraints on the tracking signal y:
%
%   y_min <= y <= y_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - y_min (optional): ny column vector, lower bound constraint values on 
%   the tracking signal
%   - y_max (optional): ny column vector, upper bound constraint values on 
%   the tracking signal
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_output_cnstr(mpc,y_min,y_max)

y_cnstr.min = y_min;
y_cnstr.max = y_max;

% Outputs box constraints
[y_cnstr.grad_min,y_cnstr.grad_max] = genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
    mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

if ~isempty(y_cnstr.min)
    y_cnstr.fi_min_x0 = zeros(mpc.Ny,1);
    [y_cnstr.hess_min,mi] = genHessIneq(y_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    y_cnstr.grad_min_feas_slv = [y_cnstr.grad_min;-ones(1,mpc.Ny)];
    
    y_cnstr.hess_min_feas_slv = genHessIneq(y_cnstr.grad_min_feas_slv);

end

if ~isempty(y_cnstr.max)
    y_cnstr.fi_max_x0 = zeros(mpc.Ny,1);
    [y_cnstr.hess_max,mi] = genHessIneq(y_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    y_cnstr.grad_max_feas_slv = [y_cnstr.grad_max;-ones(1,mpc.Ny)];
    
    y_cnstr.hess_max_feas_slv = genHessIneq(y_cnstr.grad_max_feas_slv);

end

mpc.y_cnstr = y_cnstr;
end