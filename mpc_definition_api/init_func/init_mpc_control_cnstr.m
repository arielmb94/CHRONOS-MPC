% INIT_MPC_CONTROL_CNSTR Add hard bounds on the control action.
%
%   mpc = INIT_MPC_CONTROL_CNSTR(mpc, u_min, u_max) constrains the control
%   actions to u_min <= u_k <= u_max. Control-action bounds are hard; no
%   slack is allowed. Use [] for a bound that is not needed. Bounds may be
%   scalars, nu-by-1 vectors, or time-varying nu-by-L matrices, where L is
%   the number of supplied horizon stages. A scalar is applied to every
%   control input and stage. If L < N, the last supplied column is reused
%   for the remaining stages.
%
%   Call this function after INIT_MPC_DYNAMICS has defined nu and before
%   calling BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc    - CHRONOS MPC structure.
%     u_min  - Lower bound: scalar, nu-by-1, or nu-by-L. Use [] for no
%              lower bound.
%     u_max  - Upper bound: scalar, nu-by-1, or nu-by-L. Use [] for no
%              upper bound.
%
%   Output:
%     mpc    - Updated CHRONOS MPC structure.
%
%   Example - limit every control input to the interval [-1, 1]:
%
%       mpc = init_mpc_control_cnstr(mpc, -1, 1);
function mpc = init_mpc_control_cnstr(mpc,u_min,u_max)
arguments
    mpc
    u_min = [];
    u_max = [];
end

if isempty(u_min) && isempty(u_max)
    return;
end

mpc.has_u_cnstr = 1;

% INPUT DIMENSION VALIDATION 
validate_column_vector(u_min, mpc.nu, 'u_min');
validate_column_vector(u_max, mpc.nu, 'u_max');

% Expand scalars to full local vectors if needed
if isscalar(u_min), u_min = u_min * ones(mpc.nu, 1); end
if isscalar(u_max), u_max = u_max * ones(mpc.nu, 1); end

u_cnstr.use_k0 = [];
u_cnstr.use_ter = [];
u_cnstr.rows_k0 = [];
u_cnstr.rows_ter = [];
u_cnstr.min_ineqRow_0 = [];
u_cnstr.min_ineqRow_k = [];
u_cnstr.max_ineqRow_0 = [];
u_cnstr.max_ineqRow_k = [];

if ~isempty(u_min)

    u_cnstr.min_limit = 1;

    mpc.ng_k(1:2) = mpc.ng_k(1:2) + mpc.nu;

    u_cnstr.min = zeros(mpc.nu,mpc.N);
    u_cnstr.min = fill_vec(u_cnstr.min, u_min, 1);

else
    u_cnstr.min_limit = [];
    u_cnstr.min = [];
end

if ~isempty(u_max)

    u_cnstr.max_limit = 1;

    mpc.ng_k(1:2) = mpc.ng_k(1:2) + mpc.nu;
    

    u_cnstr.max = zeros(mpc.nu,mpc.N);
    u_cnstr.max = fill_vec(u_cnstr.max, u_max, 1);
else
    u_cnstr.max_limit = [];
    u_cnstr.max = [];
end

mpc.u_cnstr = u_cnstr;

end
