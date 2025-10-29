%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_state_cnstr(mpc,x_min,x_max)
%
% Define constraints on the state vector x:
%
%   x_min <= x <= x_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - x_min (optional): nx column vector, lower bound constraint values on 
%   the state vector
%   - x_max (optional): nx column vector, upper bound constraint values on 
%   the state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_state_cnstr(mpc,x_min,x_max,x_min_activ,x_max_activ)

s_cnstr.min = x_min;
s_cnstr.max = x_max;

[s_cnstr.grad_min,s_cnstr.grad_max] = genGradX(mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

if ~isempty(s_cnstr.min)
    s_cnstr.fi_min_x0 = zeros(mpc.Nx,1);
    [s_cnstr.hess_min,mi] = genHessIneq(s_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    s_cnstr.grad_min_feas_slv = [s_cnstr.grad_min;-ones(1,mpc.Nx)];
    
    s_cnstr.hess_min_feas_slv = genHessIneq(s_cnstr.grad_min_feas_slv);

end

if ~isempty(s_cnstr.max)
    s_cnstr.fi_max_x0 = zeros(mpc.Nx,1);
    [s_cnstr.hess_max,mi] = genHessIneq(s_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    s_cnstr.grad_max_feas_slv = [s_cnstr.grad_max;-ones(1,mpc.Nx)];
    
    s_cnstr.hess_max_feas_slv = genHessIneq(s_cnstr.grad_max_feas_slv);
    
end

% Active-set like optimization

if exist("x_min_activ")

    if ~isempty(x_min_activ)

        % check that activation range is within variable range
        index = x_min_activ < x_min;
        x_min_activ(index) = x_min(index);

        s_cnstr.min_activ_set = 1;
        s_cnstr.min_activ_lim = x_min_activ;
        s_cnstr.min_activ_indicator = zeros(mpc.nx,1);
    else
        s_cnstr.min_activ_set = 0;
        s_cnstr.min_activ_lim = [];
        s_cnstr.min_activ_indicator = [];
    end
else
    s_cnstr.min_activ_set = 0;
    s_cnstr.min_activ_lim = [];
    s_cnstr.min_activ_indicator = [];
end

if exist("x_max_activ")

    if ~isempty(x_max_activ)

        % check that activation range is within variable range
        index = x_max_activ > x_max;
        x_max_activ(index) = x_max(index);

        s_cnstr.max_activ_set = 1;
        s_cnstr.max_activ_lim = x_max_activ;
        s_cnstr.max_activ_indicator = zeros(mpc.nx,1);
    else
        s_cnstr.max_activ_set = 0;
        s_cnstr.max_activ_lim = [];
        s_cnstr.max_activ_indicator = [];
    end
else
    s_cnstr.max_activ_set = 0;
    s_cnstr.max_activ_lim = [];
    s_cnstr.max_activ_indicator = [];
end

mpc.s_cnstr = s_cnstr;
end