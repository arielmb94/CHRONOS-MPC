% INIT_MPC_CONTROL_RATE_CNSTR Add bounds on the control-action rate.
%
%   mpc = INIT_MPC_CONTROL_RATE_CNSTR(mpc, du_min, du_max) constrains
%   du_min <= Delta_u_k <= du_max, where Delta_u_0 = u_0 - u_prev and
%   Delta_u_k = u_k - u_(k-1). Control-rate bounds are hard and do not
%   allow slacks. Use [] for a bound that is not needed. Bounds may be
%   scalars, nu-by-1 vectors, or time-varying nu-by-L matrices, where L is
%   the number of supplied horizon stages. A scalar is applied to every
%   control input and stage. If L < N, the last supplied column is reused
%   for the remaining stages.
%
%   Call this function after INIT_MPC_DYNAMICS has defined nu and before
%   calling BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     du_min  - Lower rate bound: scalar, nu-by-1, or nu-by-L. Use [] for
%               no lower bound.
%     du_max  - Upper rate bound: scalar, nu-by-1, or nu-by-L. Use [] for
%               no upper bound.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - limit each control move to +/-0.1 per sample:
%
%       mpc = init_mpc_control_rate_cnstr(mpc, -0.1, 0.1);
function mpc = init_mpc_control_rate_cnstr(mpc,du_min,du_max)
arguments
    mpc
    du_min = [];
    du_max = [];
end

if isempty(du_min) && isempty(du_max)
    return;
end

mpc.has_du_cnstr = 1;
mpc.has_du = 1;

% INPUT DIMENSION VALIDATION 
validate_column_vector(du_min, mpc.nu, 'du_min');
validate_column_vector(du_max, mpc.nu, 'du_max');

% Expand scalars to full local vectors if needed
if isscalar(du_min), du_min = du_min * ones(mpc.nu, 1); end
if isscalar(du_max), du_max = du_max * ones(mpc.nu, 1); end

du_cnstr.use_k0 = [];
du_cnstr.use_ter = [];
du_cnstr.rows_k0 = [];
du_cnstr.rows_ter = [];
du_cnstr.min_ineqRow_0 = [];
du_cnstr.min_ineqRow_k = [];
du_cnstr.max_ineqRow_0 = [];
du_cnstr.max_ineqRow_k = [];

if ~isempty(du_min)

    du_cnstr.min_limit = 1;
    
    mpc.ng_k(1:2) = mpc.ng_k(1:2) + mpc.nu;

    du_cnstr.min = zeros(mpc.nu,mpc.N);
    du_cnstr.min = fill_vec(du_cnstr.min, du_min, 1);
else
    du_cnstr.min_limit = [];
    du_cnstr.min = [];
end

if ~isempty(du_max)

    du_cnstr.max_limit = 1;

    mpc.ng_k(1:2) = mpc.ng_k(1:2) + mpc.nu;


    du_cnstr.max = zeros(mpc.nu,mpc.N);
    du_cnstr.max = fill_vec(du_cnstr.max, du_max, 1);    
else
    du_cnstr.max_limit = [];
    du_cnstr.max = [];
end

mpc.du_cnstr = du_cnstr;

end
