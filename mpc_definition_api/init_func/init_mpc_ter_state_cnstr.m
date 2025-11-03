%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max)
%
% Define box constraints only on the final state vector xN of the MPC 
% prediction horizon:
%
%   xN_min <= xN <= xN_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - x_ter_min (optional): nx column vector, lower bound constraint values
%   on the terminal state vector
%   - x_ter_max (optional): nx column vector, upper bound constraint values
%   on the terminal state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max,...
                                            x_ter_min_activ,x_ter_max_activ)

s_ter_cnstr.min = x_ter_min; 
s_ter_cnstr.max = x_ter_max;

[s_ter_cnstr.grad_min,s_ter_cnstr.grad_max] = genGradXter(mpc.N,mpc.N_ctr_hor,...
    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(s_ter_cnstr.min)
    s_ter_cnstr.fi_min_x0 = zeros(mpc.nx,1);
    [s_ter_cnstr.hess_min,mi] = genHessIneq(s_ter_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    s_ter_cnstr.grad_min_feas_slv = [s_ter_cnstr.grad_min;-ones(1,mpc.nx)];
    
    s_ter_cnstr.hess_min_feas_slv = genHessIneq(s_ter_cnstr.grad_min_feas_slv);

    % Active-set like optimization setting
    if exist("x_ter_min_activ")
        s_ter_cnstr = set_active_set_min(s_ter_cnstr,x_ter_min,x_ter_min_activ,mpc.nx);
    else
        s_ter_cnstr = empty_active_set_min(s_ter_cnstr);
    end
end

if ~isempty(s_ter_cnstr.max)
    s_ter_cnstr.fi_max_x0 = zeros(mpc.nx,1);
    [s_ter_cnstr.hess_max,mi] = genHessIneq(s_ter_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    s_ter_cnstr.grad_max_feas_slv = [s_ter_cnstr.grad_max;-ones(1,mpc.nx)];
    
    s_ter_cnstr.hess_max_feas_slv = genHessIneq(s_ter_cnstr.grad_max_feas_slv);

    % Active-set like optimization setting
    if exist("x_ter_max_activ")
        s_ter_cnstr = set_active_set_max(s_ter_cnstr,x_ter_max,x_ter_max_activ,mpc.nx);
    else
        s_ter_cnstr = empty_active_set_max(s_ter_cnstr);
    end
end

mpc.s_ter_cnstr = s_ter_cnstr;
end