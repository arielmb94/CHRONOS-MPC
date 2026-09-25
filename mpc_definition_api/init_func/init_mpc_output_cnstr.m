% INIT_MPC_OUTPUT_CNSTR Add bounds on the predicted output.
%
%   mpc = INIT_MPC_OUTPUT_CNSTR(mpc, y_min, y_max) constrains the output
%   signal defined by INIT_MPC_DYNAMICS (default: y = s) or
%   INIT_MPC_OUTPUT. Use [] for a bound that is not needed. Bounds may be
%   scalars, ny-by-1 vectors, or time-varying ny-by-L matrices, where L is
%   the number of supplied horizon stages. A scalar is applied to every
%   output and stage. If L < N, the last supplied column is reused for the
%   remaining stages.
%
%   Output bounds apply at interior stages k = 1,...,N-1. They also apply at
%   k = 0 to output rows that depend on u, and at k = N to output rows that
%   depend exclusively on the state, as determined when the output model is
%   initialized.
%
%   mpc = INIT_MPC_OUTPUT_CNSTR(..., qv_min, qv_max) also sets the linear
%   penalties for lower- and upper-bound violations. Output constraints are
%   soft: CHRONOS may violate a bound through a feasibility slack when the
%   bound cannot be satisfied. Larger qv values make violations more costly.
%   Penalties use the same scalar, vector, or time-varying horizon layout as
%   the bounds. Leave a penalty empty to let CHRONOS select its default
%   during BUILD_CHRONOS_MPC.
%
%   Call this function after defining the output model and before calling
%   BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     y_min   - Lower bound: scalar, ny-by-1, or ny-by-L. Use [] for no
%               lower bound.
%     y_max   - Upper bound: scalar, ny-by-1, or ny-by-L. Use [] for no
%               upper bound.
%     qv_min  - Optional lower-bound violation penalty: scalar, ny-by-1,
%               or time-varying ny-by-L.
%     qv_max  - Optional upper-bound violation penalty: scalar, ny-by-1,
%               or time-varying ny-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - constrain a scalar output with soft bounds:
%
%       C = [1, 0];
%       D = zeros(1, mpc.nu);
%       mpc = init_mpc_output(mpc, C, D);
%       mpc = init_mpc_output_cnstr(mpc, -2, 2, 100, 100);
function mpc = init_mpc_output_cnstr(mpc,y_min,y_max,qv_min,qv_max)
arguments
    mpc
    y_min = []
    y_max = []
    qv_min = []
    qv_max = []
end

if isempty(y_min) && isempty(y_max)
    return;
end

% A penalty is ignored when its bound is inactive or when the supplied
% penalty is all zero; build_chronos_mpc will select the default penalty.
if isempty(y_min) || isempty(qv_min) || ~any(qv_min(:)), qv_min = []; end
if isempty(y_max) || isempty(qv_max) || ~any(qv_max(:)), qv_max = []; end

mpc.has_y_cnstr = 1;

% INPUT DIMENSION VALIDATION 
validate_column_vector(y_min, mpc.ny, 'y_min');
validate_column_vector(y_max, mpc.ny, 'y_max');
validate_column_vector(qv_min, mpc.ny, 'qv_min');
validate_column_vector(qv_max, mpc.ny, 'qv_max');

y_cnstr.use_k0 = mpc.y_use_k0;
y_cnstr.use_ter = mpc.y_use_ter;
y_cnstr.rows_k0 = mpc.y_rows_k0;
y_cnstr.rows_ter = mpc.y_rows_ter;
y_cnstr.use_s = mpc.y_use_s;
y_cnstr.use_u = mpc.y_use_u;
y_cnstr.use_d = mpc.y_use_d;

y_cnstr.min_ineqRow_0 = [];
y_cnstr.min_ineqRow_k = [];
y_cnstr.min_ineqRow_ter = [];
y_cnstr.max_ineqRow_0 = [];
y_cnstr.max_ineqRow_k = [];
y_cnstr.max_ineqRow_ter = [];
y_cnstr.min_row_v_0 = [];
y_cnstr.min_row_v_k = [];
y_cnstr.max_row_v_0 = [];
y_cnstr.max_row_v_k = [];

% Expand scalars to full local vectors if needed
if isscalar(y_min), y_min = y_min * ones(mpc.ny, 1); end
if isscalar(y_max), y_max = y_max * ones(mpc.ny, 1); end
if isscalar(qv_min), qv_min = qv_min * ones(mpc.ny, 1); end
if isscalar(qv_max), qv_max = qv_max * ones(mpc.ny, 1); end

if ~isempty(y_min)

    y_cnstr.min_limit = 1;
    
    if ~isempty(mpc.y_use_k0)
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.ny_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.ny_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.ny;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.ny;
    if ~isempty(mpc.y_use_ter)
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.ny_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.ny_ter;
    end

    y_full = zeros(mpc.ny, mpc.N);
    y_full = fill_vec(y_full, y_min, 1);
    y_cnstr.min = y_full(:,1:mpc.N-1);
    if ~isempty(mpc.y_use_k0)
        y_cnstr.min_0 = y_full(mpc.y_rows_k0,1);
    else
        y_cnstr.min_0 = [];
    end
    if ~isempty(mpc.y_use_ter)
        y_cnstr.min_ter = y_full(mpc.y_rows_ter,mpc.N);
    else
        y_cnstr.min_ter = [];
    end

    qv_min_full = zeros(mpc.ny, mpc.N);
    if ~isempty(qv_min)
        qv_min_full = fill_vec(qv_min_full, qv_min, 1);
    end
    y_cnstr.qv_min = qv_min_full(:,1:mpc.N-1);
    if ~isempty(mpc.y_use_k0)
        y_cnstr.qv_min_0 = qv_min_full(mpc.y_rows_k0,1);
    else
        y_cnstr.qv_min_0 = [];
    end
    if ~isempty(mpc.y_use_ter)
        y_cnstr.qv_min_ter = qv_min_full(mpc.y_rows_ter,mpc.N);
    else
        y_cnstr.qv_min_ter = [];
    end

else
    y_cnstr.min_limit = [];
    y_cnstr.min = [];
    y_cnstr.min_0 = [];
    y_cnstr.min_ter = [];
    y_cnstr.qv_min = [];
    y_cnstr.qv_min_0 = [];
    y_cnstr.qv_min_ter = [];
end

if ~isempty(y_max)

    y_cnstr.max_limit = 1;
    
    if ~isempty(mpc.y_use_k0)
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.ny_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.ny_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.ny;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.ny;
    if ~isempty(mpc.y_use_ter)
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.ny_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.ny_ter;
    end

    y_full = zeros(mpc.ny, mpc.N);
    y_full = fill_vec(y_full, y_max, 1);
    y_cnstr.max = y_full(:,1:mpc.N-1);
    if ~isempty(mpc.y_use_k0)
        y_cnstr.max_0 = y_full(mpc.y_rows_k0,1);
    else
        y_cnstr.max_0 = [];
    end
    if ~isempty(mpc.y_use_ter)
        y_cnstr.max_ter = y_full(mpc.y_rows_ter,mpc.N);
    else
        y_cnstr.max_ter = [];
    end

    qv_max_full = zeros(mpc.ny, mpc.N);
    if ~isempty(qv_max)
        qv_max_full = fill_vec(qv_max_full, qv_max, 1);
    end
    y_cnstr.qv_max = qv_max_full(:,1:mpc.N-1);
    if ~isempty(mpc.y_use_k0)
        y_cnstr.qv_max_0 = qv_max_full(mpc.y_rows_k0,1);
    else
        y_cnstr.qv_max_0 = [];
    end
    if ~isempty(mpc.y_use_ter)
        y_cnstr.qv_max_ter = qv_max_full(mpc.y_rows_ter,mpc.N);
    else
        y_cnstr.qv_max_ter = [];
    end

else
    y_cnstr.max_limit = [];
    y_cnstr.max = [];
    y_cnstr.max_0 = [];
    y_cnstr.max_ter = [];
    y_cnstr.qv_max = [];
    y_cnstr.qv_max_0 = [];
    y_cnstr.qv_max_ter = [];
end

mpc.y_cnstr = y_cnstr;

end
