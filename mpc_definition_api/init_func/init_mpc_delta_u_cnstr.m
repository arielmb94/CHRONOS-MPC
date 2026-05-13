% INIT_MPC_DELTA_U_CNSTR Defines control action rate constraints.
%
%   mpc = INIT_MPC_DELTA_U_CNSTR(mpc, du_min, du_max) sets strict (hard) lower and 
%   upper bounds on the control action rate. The solver will strictly enforce these 
%   limits.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       du_min              - [nu x 1] Array of lower control action rate limits (use [] if none).
%       du_max              - [nu x 1] Array of upper control action rate limits (use [] if none).
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (constraint gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
function mpc = init_mpc_delta_u_cnstr(mpc,du_min,du_max)
arguments
    mpc
    du_min = [];
    du_max = [];
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(du_min, mpc.nu, 'du_min');
validate_column_vector(du_max, mpc.nu, 'du_max');

% Expand scalars to full local vectors if needed
if isscalar(du_min), du_min = du_min * ones(mpc.nu, 1); end
if isscalar(du_max), du_max = du_max * ones(mpc.nu, 1); end

du_cnstr.min = du_min;
du_cnstr.max = du_max;

if ~isempty(du_cnstr.min)

    du_cnstr.min_limit = 1;

    du_cnstr.fi_min_x0 = zeros(mpc.Nu,1);

    % Differential Control box constraints
    du_cnstr.grad_min = -1 * genGradDeltaU(mpc.N_ctr_hor,...
                             mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    [du_cnstr.hess_min,mi] = genHessIneq(du_cnstr.grad_min);
    mpc.m = mpc.m+mi;
else
    du_cnstr.min_limit = 0;
end

if ~isempty(du_cnstr.max)

    du_cnstr.max_limit = 1;
    
    du_cnstr.fi_max_x0 = zeros(mpc.Nu,1);

    % Differential Control box constraints
    du_cnstr.grad_max = genGradDeltaU(mpc.N_ctr_hor,...
                        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    [du_cnstr.hess_max,mi] = genHessIneq(du_cnstr.grad_max);
    mpc.m = mpc.m+mi;
else
    du_cnstr.max_limit = 0;
end

mpc.du_cnstr = du_cnstr;

end