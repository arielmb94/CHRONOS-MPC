%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   %mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,h_min,h_max)
%
% Define constraints on a custom user defined signal h. The signal h is
% defined as:
%
%   h = Ch * x + Dh * u + Ddh * dh
%
% The user must define the signal h by selecting appropiate values for the
% matrices Ch, Dh, Ddh.
% 
% The constraints on h are then of the form:
%
%   h_min <= h <= h_max
%
% In:
%   - mpc: CHRONOS mpc structure
%   - h_min: nh column vector, lower bound constraint values on the user
%   defined signal h
%   - h_max: nh column vector, upper bound constraint values on the user
%   defined signal h
%   - Ch: nh x nx matrix, states output matrix
%   - Dh: nh x nu matrix, input feedtrhough matrix
%   - Ddh: nh x ndh matrix, disturbance feedtrhough matrix
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
% Example:
% We want to limit the variation of the MPC control action with respect a
% a given value u_star, e.g. we want to constraint h = u - u_star.
% To achieve this we define h by selecting:
%   - Ch = [0 0 ... 0]
%   - Dh = [1]
%   - Ddh = [-1]
% which corresponds to h = [0 0 ... 0] * x + [1] * u + [-1] * dh
% The "disturbance" term on h corresponds to u_star, to be introduced on 
% the appropiate field on mpc_solve() during runtime MPC execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Ddh,...
                                        h_min,h_max,h_min_activ,h_max_activ)

h_cnstr.min = h_min;
h_cnstr.max = h_max;

% General Inequality Matrix
mpc.Ch = Ch;
mpc.Dh = Dh;
mpc.Ddh = Ddh;

mpc.ndh = size(Ddh,2);  %number of disturbance inputs to general inequalities
mpc.nh = size(Ch,1);  %number of general inequalities

if mpc.Dh == 0
    mpc.Nh = (mpc.N-1)*mpc.nh;
else
    mpc.Nh = mpc.N*mpc.nh;
end

mpc.h = zeros(mpc.Nh,1);

mpc.Ndh = mpc.N*mpc.ndh;

% General Inequalites box constraints
[h_cnstr.grad_min,h_cnstr.grad_max] = genGradY(mpc.Ch,mpc.Dh,mpc.N,mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.Nh,mpc.nx,mpc.nu,mpc.nh);

if ~isempty(h_cnstr.min)
    h_cnstr.fi_min_x0 = zeros(mpc.Nh,1);
    [h_cnstr.hess_min,mi] = genHessIneq(h_cnstr.grad_min);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver min condition
    h_cnstr.grad_min_feas_slv = [h_cnstr.grad_min;-ones(1,mpc.Nh)];
    
    h_cnstr.hess_min_feas_slv = genHessIneq(h_cnstr.grad_min_feas_slv);

end
if ~isempty(h_cnstr.max)
    h_cnstr.fi_max_x0 = zeros(mpc.Nh,1);
    [h_cnstr.hess_max,mi] = genHessIneq(h_cnstr.grad_max);
    mpc.m = mpc.m+mi;

    % initialize feasibility solver max condition
    h_cnstr.grad_max_feas_slv = [h_cnstr.grad_max;-ones(1,mpc.Nh)];
    
    h_cnstr.hess_max_feas_slv = genHessIneq(h_cnstr.grad_max_feas_slv);

end

% Active-set like optimization

if exist("h_min_activ")

    if ~isempty(h_min_activ)

        % check that activation range is within variable range
        index = h_min_activ < h_min;
        h_min_activ(index) = h_min(index);

        h_cnstr.min_activ_set = 1;
        h_cnstr.min_activ_lim = h_min_activ;
        h_cnstr.min_activ_indicator = zeros(mpc.nh,1);
    else
        h_cnstr.min_activ_set = 0;
        h_cnstr.min_activ_lim = [];
        h_cnstr.min_activ_indicator = [];
    end
else
    h_cnstr.min_activ_set = 0;
    h_cnstr.min_activ_lim = [];
    h_cnstr.min_activ_indicator = [];
end

if exist("h_max_activ")

    if ~isempty(h_max_activ)

        % check that activation range is within variable range
        index = h_max_activ > h_max;
        h_max_activ(index) = h_max(index);

        h_cnstr.max_activ_set = 1;
        h_cnstr.max_activ_lim = h_max_activ;
        h_cnstr.max_activ_indicator = zeros(mpc.nh,1);
    else
        h_cnstr.max_activ_set = 0;
        h_cnstr.max_activ_lim = [];
        h_cnstr.max_activ_indicator = [];
    end
else
    h_cnstr.max_activ_set = 0;
    h_cnstr.max_activ_lim = [];
    h_cnstr.max_activ_indicator = [];
end

mpc.h_cnstr = h_cnstr;

end