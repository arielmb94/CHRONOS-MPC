function mpc = init_mpc_state_cnstr(mpc,x_min,x_max,qv_min,qv_max)
% INIT_MPC_STATE_CNSTR Add bounds on the predicted state.
%
%   mpc = INIT_MPC_STATE_CNSTR(mpc, x_min, x_max) constrains the predicted
%   states to x_min <= x_k <= x_max. Use [] for a bound that is not needed.
%   Bounds may be scalars, nx-by-1 vectors, or time-varying nx-by-L
%   matrices, where L is the number of supplied horizon stages. A scalar is
%   applied to every state and stage. If L < N, the last supplied column is
%   reused for the remaining stages.
%
%   mpc = INIT_MPC_STATE_CNSTR(..., qv_min, qv_max) also sets the linear
%   penalties for lower- and upper-bound violations. State constraints are
%   soft: CHRONOS may violate a bound through a feasibility slack when the
%   bound cannot be satisfied. Larger qv values make violations more costly.
%   Penalties use the same scalar, vector, or time-varying horizon layout as
%   the bounds. Leave a penalty empty to let CHRONOS select its default
%   during BUILD_CHRONOS_MPC.
%
%   Call this function after INIT_MPC_DYNAMICS and before calling 
%   BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     x_min   - Lower bound: scalar, nx-by-1, or nx-by-L. Use [] for no
%               lower bound.
%     x_max   - Upper bound: scalar, nx-by-1, or nx-by-L. Use [] for no
%               upper bound.
%     qv_min  - Optional lower-bound violation penalty: scalar, nx-by-1,
%               or time-varying nx-by-L.
%     qv_max  - Optional upper-bound violation penalty: scalar, nx-by-1,
%               or time-varying nx-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - apply bounds to a two-state model:
%
%       mpc = init_mpc_state_cnstr(mpc, [-5; -2], [5; 2], 100, 100);

arguments
    mpc
    x_min = []
    x_max = []
    qv_min = []
    qv_max = []
end

if isempty(x_min) && isempty(x_max)
    return;
end

% A penalty is ignored when its bound is inactive or when the supplied
% penalty is all zero; build_chronos_mpc will select the default penalty.
if isempty(x_min) || isempty(qv_min) || ~any(qv_min(:)), qv_min = []; end
if isempty(x_max) || isempty(qv_max) || ~any(qv_max(:)), qv_max = []; end

mpc.has_x_cnstr = 1;

% INPUT DIMENSION VALIDATION 
validate_column_vector(x_min, mpc.nx, 'x_min');
validate_column_vector(x_max, mpc.nx, 'x_max');
validate_column_vector(qv_min, mpc.nx, 'qv_min');
validate_column_vector(qv_max, mpc.nx, 'qv_max');

% Expand scalars to full local vectors if needed
if isscalar(x_min), x_min = x_min * ones(mpc.nx, 1); end
if isscalar(x_max), x_max = x_max * ones(mpc.nx, 1); end
if isscalar(qv_min), qv_min = qv_min * ones(mpc.nx, 1); end
if isscalar(qv_max), qv_max = qv_max * ones(mpc.nx, 1); end

x_cnstr.use_k0 = [];
x_cnstr.use_ter = 1;
x_cnstr.rows_k0 = [];
x_cnstr.rows_ter = 1:mpc.nx;
x_cnstr.min_ineqRow_k = [];
x_cnstr.min_ineqRow_ter = [];
x_cnstr.max_ineqRow_k = [];
x_cnstr.max_ineqRow_ter = [];
x_cnstr.min_row_v_k = [];
x_cnstr.max_row_v_k = [];

if ~isempty(x_min)

    x_cnstr.min_limit = 1;
   
    mpc.ng_k(2:3) = mpc.ng_k(2:3) + mpc.nx;
    mpc.nv_k(2:3) = mpc.nv_k(2:3) + mpc.nx;

    x_cnstr.min = zeros(mpc.nx,mpc.N);
    x_cnstr.min = fill_vec(x_cnstr.min, x_min, 1);
    
    % Initialize Penalty term for new slack variables
    % Allocate and fill time-varying slack penalties
    x_cnstr.qv_min = zeros(mpc.nx, mpc.N-1);
    x_cnstr.qv_min_ter = zeros(mpc.nx, 1);
    if ~isempty(qv_min)
        x_cnstr.qv_min = fill_vec(x_cnstr.qv_min, qv_min, 1);
        x_cnstr.qv_min_ter = qv_min(:,min(size(qv_min,2),mpc.N));
    end
else
    x_cnstr.min_limit = [];
    x_cnstr.min = [];
    x_cnstr.qv_min = [];
    x_cnstr.qv_min_ter = [];
end

if ~isempty(x_max)

    x_cnstr.max_limit = 1;
   
    mpc.ng_k(2:3) = mpc.ng_k(2:3) + mpc.nx;
    mpc.nv_k(2:3) = mpc.nv_k(2:3) + mpc.nx;

    x_cnstr.max = zeros(mpc.nx,mpc.N);
    x_cnstr.max = fill_vec(x_cnstr.max, x_max, 1);
    
    % Initialize Penalty term for new slack variables
    % Allocate and fill time-varying slack penalties
    x_cnstr.qv_max = zeros(mpc.nx, mpc.N-1);
    x_cnstr.qv_max_ter = zeros(mpc.nx, 1);
    if ~isempty(qv_max)
        x_cnstr.qv_max = fill_vec(x_cnstr.qv_max, qv_max, 1);
        x_cnstr.qv_max_ter = qv_max(:,min(size(qv_max,2),mpc.N));
    end
else
    x_cnstr.max_limit = [];
    x_cnstr.max = [];
    x_cnstr.qv_max = [];
    x_cnstr.qv_max_ter = [];
end

mpc.x_cnstr = x_cnstr;

end
