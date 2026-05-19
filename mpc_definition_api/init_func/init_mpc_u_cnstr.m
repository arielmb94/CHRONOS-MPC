% INIT_MPC_U_CNSTR Defines control action constraints and soft-constraint penalties.
%
%   mpc = INIT_MPC_U_CNSTR(mpc, u_min, u_max) sets strict (hard) lower and 
%   upper bounds on the control action. The solver will strictly enforce these 
%   limits.
%
%   mpc = INIT_MPC_U_CNSTR(mpc, u_min, u_max, u_min_slack_active, u_max_slack_active, qv_min, qv_max) 
%   allows you to define specific bounds as "soft" constraints. For control
%   actions it is advisable to use hard constraints since it is a magnitude
%   the solver will have direct control over.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure
%       u_min              - [nu x 1] Array of lower control action limits (use [] if none).
%       u_max              - [nu x 1] Array of upper control action limits (use [] if none).
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (constraint gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
function mpc = init_mpc_u_cnstr(mpc,u_min,u_max)
arguments
    mpc
    u_min = [];
    u_max = [];
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(u_min, mpc.nu, 'u_min');
validate_column_vector(u_max, mpc.nu, 'u_max');

% Expand scalars to full local vectors if needed
if isscalar(u_min), u_min = u_min * ones(mpc.nu, 1); end
if isscalar(u_max), u_max = u_max * ones(mpc.nu, 1); end

u_cnstr.min = u_min;
u_cnstr.max = u_max;

if ~isempty(u_cnstr.min)

    u_cnstr.min_limit = 1;

    u_cnstr.fi_min_x0 = zeros(mpc.Nu,1);

    % Control box constraints
    u_cnstr.grad_min = -1 * genGradU(mpc.N_ctr_hor,...
                                     mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    [u_cnstr.hess_min,mi] = genHessIneq(u_cnstr.grad_min);
    mpc.m = mpc.m+mi;
else
    u_cnstr.min_limit = 0;
end

if ~isempty(u_cnstr.max)

    u_cnstr.max_limit = 1;

    u_cnstr.fi_max_x0 = zeros(mpc.Nu,1);

    % Control box constraints
    u_cnstr.grad_max = genGradU(mpc.N_ctr_hor,...
                                mpc.Nx,mpc.Nu,mpc.nx,mpc.nu,mpc.Nv);

    [u_cnstr.hess_max,mi] = genHessIneq(u_cnstr.grad_max);
    mpc.m = mpc.m+mi;
    
else
    u_cnstr.max_limit = 0;
end

mpc.u_cnstr = u_cnstr;

end