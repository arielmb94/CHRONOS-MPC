% INIT_MPC_CUSTOM_CNSTR Add bounds on a user-defined linear signal.
%
%   mpc = INIT_MPC_CUSTOM_CNSTR(mpc, Ch, Dh, Dsuh, Ddh, h_min, h_max)
%   defines the custom signal
%
%       h_k = Ch_k*s_k + Dh_k*u_k + Dsuh_k*su_k + Ddh_k*dh_k
%
%   and constrains it to h_min <= h_k <= h_max. Here, su_k is the control
%   action preceding u_k, and dh_k is a dedicated fixed known input for the
%   custom constraint. Use [] for a coefficient or bound that is not needed.
%   Coefficient matrices, bounds, and penalties may vary over the horizon.
%   L denotes the number of supplied horizon stages; if L < N, the last
%   supplied stage is reused for the remaining stages.
%
%   mpc = INIT_MPC_CUSTOM_CNSTR(..., qv_min, qv_max) also sets the linear
%   penalties for lower- and upper-bound violations. Custom constraints are
%   soft: CHRONOS may violate a bound through a feasibility slack when the
%   bound cannot be satisfied. Larger qv values make violations more costly.
%   Leave a penalty empty to let CHRONOS select its default during
%   BUILD_CHRONOS_MPC.
%
%   Call this function after defining the MPC model dimensions and before
%   calling BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     Ch      - State coefficient, size nh-by-nx or nh-by-nx-by-L.
%     Dh      - Control coefficient, size nh-by-nu or nh-by-nu-by-L.
%     Dsuh    - Previous-control coefficient, size nh-by-nu or nh-by-nu-by-L.
%     Ddh     - Fixed-known-input coefficient, size nh-by-ndh or
%               nh-by-ndh-by-L.
%     h_min   - Lower bound: scalar, nh-by-1, or nh-by-L. Use [] for no
%               lower bound.
%     h_max   - Upper bound: scalar, nh-by-1, or nh-by-L. Use [] for no
%               upper bound.
%     qv_min  - Optional lower-bound violation penalty: scalar, nh-by-1,
%               or time-varying nh-by-L.
%     qv_max  - Optional upper-bound violation penalty: scalar, nh-by-1,
%               or time-varying nh-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - constrain the deviation from a runtime control reference:
%
%       % h = u - u_ref, with -0.2 <= h <= 0.2
%       Ch   = zeros(nu, nx);
%       Dh   = eye(nu);
%       Dsuh = [];
%       Ddh  = -eye(nu);
%       mpc = init_mpc_custom_cnstr(mpc, Ch, Dh, Dsuh, Ddh, ...
%                                    -0.2, 0.2, 100, 100);
%
%   Pass u_ref as the dh_in argument of MPC_SOLVE at runtime.
function mpc = init_mpc_custom_cnstr(mpc,Ch,Dh,Dsuh,Ddh,...
                                            h_min,h_max, ...
                                            qv_min,qv_max)
arguments
    mpc
    Ch = []
    Dh = []
    Dsuh = []
    Ddh = []
    h_min = []
    h_max = []
    qv_min = []
    qv_max = []
end

% Empty and all-zero optional coefficients are omitted before dimension
% inference or validation, independently of their supplied dimensions.
if ~isempty(Ch) && ~any(Ch(:)), Ch = []; end
if ~isempty(Dh) && ~any(Dh(:)), Dh = []; end
if ~isempty(Dsuh) && ~any(Dsuh(:)), Dsuh = []; end
if ~isempty(Ddh) && ~any(Ddh(:)), Ddh = []; end

% A custom constraint requires an active bound and a decision-dependent
% primary signal term. Ddh alone cannot define the signal.
if (isempty(h_min) && isempty(h_max)) || ...
        (isempty(Ch) && isempty(Dh) && isempty(Dsuh))
    return;
end

% A penalty is ignored when its bound is inactive or when the supplied
% penalty is all zero; build_chronos_mpc will select the default penalty.
if isempty(h_min) || isempty(qv_min) || ~any(qv_min(:)), qv_min = []; end
if isempty(h_max) || isempty(qv_max) || ~any(qv_max(:)), qv_max = []; end

%number of general inequalities
if ~isempty(Ch)
    mpc.nh = size(Ch,1);  
elseif ~isempty(Dh)
    mpc.nh = size(Dh,1);
else
    mpc.nh = size(Dsuh,1);
end

validate_matrix(Ch, mpc.nh, mpc.nx, 'Ch');
validate_matrix(Dh, mpc.nh, mpc.nu, 'Dh');
validate_matrix(Dsuh, mpc.nh, mpc.nu, 'Dsuh');
if ~isempty(Ddh)
    mpc.ndh = size(Ddh,2);
    validate_matrix(Ddh, mpc.nh, mpc.ndh, 'Ddh');
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(h_min, mpc.nh, 'h_min');
validate_column_vector(h_max, mpc.nh, 'h_max');
validate_column_vector(qv_min, mpc.nh, 'qv_min');
validate_column_vector(qv_max, mpc.nh, 'qv_max');

% general constraints boolean
mpc.has_h_cnstr = 1;

h_cnstr.use_s = [];
h_cnstr.use_u = [];
h_cnstr.use_su = [];
h_cnstr.use_d = [];
h_cnstr.rows_k0 = [];
h_cnstr.rows_ter = [];
h_cnstr.min_ineqRow_0 = [];
h_cnstr.min_ineqRow_k = [];
h_cnstr.min_ineqRow_ter = [];
h_cnstr.max_ineqRow_0 = [];
h_cnstr.max_ineqRow_k = [];
h_cnstr.max_ineqRow_ter = [];
h_cnstr.min_row_v_0 = [];
h_cnstr.min_row_v_k = [];
h_cnstr.max_row_v_0 = [];
h_cnstr.max_row_v_k = [];

mpc.Ch_0 = [];
mpc.Dh_0 = [];
mpc.Dsuh_0 = [];
mpc.Ddh_0 = [];
mpc.Ch_ter = [];

if ~isempty(Ch)
    h_cnstr.use_s = 1;

    len_Ch = size(Ch,3);
    if len_Ch < mpc.N
        mpc.Ch = zeros(mpc.nh,mpc.nx,mpc.N-1);
        mpc.Ch = fill_mat(mpc.Ch, Ch, 1);
        Ch_ter = mpc.Ch(:,:,mpc.N-1);
    else
        mpc.Ch = Ch(:,:,1:mpc.N-1);
        Ch_ter = Ch(:,:,mpc.N);
    end
end

if ~isempty(Dh)
    h_cnstr.use_u = 1;

    len_Dh = size(Dh,3);
    if len_Dh < mpc.N
        mpc.Dh = zeros(mpc.nh,mpc.nu,mpc.N-1);
        mpc.Dh = fill_mat(mpc.Dh, Dh, 1);
        Dh_ter = mpc.Dh(:,:,mpc.N-1);
    else
        mpc.Dh = Dh(:,:,1:mpc.N-1);
        Dh_ter = Dh(:,:,mpc.N);
    end
end

if ~isempty(Dsuh)
    h_cnstr.use_su = 1;
    mpc.has_du = 1;

    len_Dsuh = size(Dsuh,3);
    if len_Dsuh < mpc.N
        mpc.Dsuh = zeros(mpc.nh,mpc.nu,mpc.N-1);
        mpc.Dsuh = fill_mat(mpc.Dsuh, Dsuh, 1);
        Dsuh_ter = mpc.Dsuh(:,:,mpc.N-1);
    else
        mpc.Dsuh = Dsuh(:,:,1:mpc.N-1);
        Dsuh_ter = Dsuh(:,:,mpc.N);
    end
end

if ~isempty(Ddh)
    h_cnstr.use_d = 1;

    len_Ddh = size(Ddh,3);
    if len_Ddh < mpc.N
        mpc.Ddh = zeros(mpc.nh,mpc.ndh,mpc.N-1);
        mpc.Ddh = fill_mat(mpc.Ddh, Ddh, 1);
        Ddh_ter = mpc.Ddh(:,:,mpc.N-1);
    else
        mpc.Ddh = Ddh(:,:,1:mpc.N-1);
        Ddh_ter = Ddh(:,:,mpc.N);
    end
end

h_cnstr.use_k0 = [];
h_cnstr.use_ter = [];
% at k = 0, only rows with Dh!=0 (with dependence on control action u) are
% considered
h_row_0 = find(~all(Dh(:,:,1)==0,2));
mpc.nh_0 = length(h_row_0);

if mpc.nh_0
    h_cnstr.rows_k0 = h_row_0;
    h_cnstr.use_k0 = 1;

    if ~isempty(h_cnstr.use_s), mpc.Ch_0 = Ch(h_cnstr.rows_k0,:,1); end
    if ~isempty(h_cnstr.use_u), mpc.Dh_0 = Dh(h_cnstr.rows_k0,:,1); end
    if ~isempty(h_cnstr.use_su), mpc.Dsuh_0 = Dsuh(h_cnstr.rows_k0,:,1); end
    if ~isempty(h_cnstr.use_d), mpc.Ddh_0 = Ddh(h_cnstr.rows_k0,:,1); end
end

% at k = N, only rows strictly dependent on s are considered
if ~isempty(h_cnstr.use_s)
    strict_s_rows = any(Ch_ter~=0,2);
    if ~isempty(h_cnstr.use_u), strict_s_rows = strict_s_rows & all(Dh_ter==0,2); end
    if ~isempty(h_cnstr.use_su), strict_s_rows = strict_s_rows & all(Dsuh_ter==0,2); end
    if ~isempty(h_cnstr.use_d), strict_s_rows = strict_s_rows & all(Ddh_ter==0,2); end

    h_row_ter = find(strict_s_rows);
    mpc.nh_ter = length(h_row_ter);
else
    mpc.nh_ter = 0;
end
if mpc.nh_ter
    h_cnstr.rows_ter = h_row_ter;
    h_cnstr.use_ter = 1;

    mpc.Ch_ter = Ch_ter(h_cnstr.rows_ter,:); 
end

% init h vector
if ~isempty(h_cnstr.use_k0), mpc.h_0 = zeros(mpc.nh_0,1); else, mpc.h_0=[]; end
mpc.h = zeros(mpc.nh,mpc.N-1);
if ~isempty(h_cnstr.use_ter), mpc.h_ter = zeros(mpc.nh_ter,1); else, mpc.h_ter=[]; end

% init dh disturbance vector
if ~isempty(h_cnstr.use_d)
    mpc.dh = zeros(mpc.ndh,mpc.N);
end

% Expand scalars to full vectors if needed
if isscalar(h_min), h_min = h_min * ones(mpc.nh, 1); end
if isscalar(h_max), h_max = h_max * ones(mpc.nh, 1); end
if isscalar(qv_min), qv_min = qv_min * ones(mpc.nh, 1); end
if isscalar(qv_max), qv_max = qv_max * ones(mpc.nh, 1); end

if ~isempty(h_min)

    h_cnstr.min_limit = 1;

    if ~isempty(h_cnstr.use_k0)
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.nh_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.nh_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.nh;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.nh;
    if ~isempty(h_cnstr.use_ter)
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.nh_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.nh_ter;
    end

    h_min_full = zeros(mpc.nh, mpc.N);
    h_min_full = fill_vec(h_min_full, h_min, 1);
    h_cnstr.min = h_min_full(:,1:mpc.N-1);
    if ~isempty(h_cnstr.use_k0)
        h_cnstr.min_0 = h_min_full(h_cnstr.rows_k0,1);
    else
        h_cnstr.min_0 = [];
    end
    if ~isempty(h_cnstr.use_ter)
        h_cnstr.min_ter = h_min_full(h_cnstr.rows_ter,mpc.N);
    else
        h_cnstr.min_ter = [];
    end

    qv_min_full = zeros(mpc.nh, mpc.N);
    if ~isempty(qv_min)
        qv_min_full = fill_vec(qv_min_full, qv_min, 1);
    end
    h_cnstr.qv_min = qv_min_full(:,1:mpc.N-1);
    if ~isempty(h_cnstr.use_k0)
        h_cnstr.qv_min_0 = qv_min_full(h_cnstr.rows_k0,1);
    else
        h_cnstr.qv_min_0 = [];
    end
    if ~isempty(h_cnstr.use_ter)
        h_cnstr.qv_min_ter = qv_min_full(h_cnstr.rows_ter,mpc.N);
    else
        h_cnstr.qv_min_ter = [];
    end
   
else
    h_cnstr.min_limit = [];
    h_cnstr.min = [];
    h_cnstr.min_0 = [];
    h_cnstr.min_ter = [];
    h_cnstr.qv_min = [];
    h_cnstr.qv_min_0 = [];
    h_cnstr.qv_min_ter = [];
end

if ~isempty(h_max)

    h_cnstr.max_limit = 1;

    if ~isempty(h_cnstr.use_k0)
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.nh_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.nh_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.nh;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.nh;
    if ~isempty(h_cnstr.use_ter)
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.nh_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.nh_ter;
    end

    h_max_full = zeros(mpc.nh, mpc.N);
    h_max_full = fill_vec(h_max_full, h_max, 1);
    h_cnstr.max = h_max_full(:,1:mpc.N-1);
    if ~isempty(h_cnstr.use_k0)
        h_cnstr.max_0 = h_max_full(h_cnstr.rows_k0,1);
    else
        h_cnstr.max_0 = [];
    end
    if ~isempty(h_cnstr.use_ter)
        h_cnstr.max_ter = h_max_full(h_cnstr.rows_ter,mpc.N);
    else
        h_cnstr.max_ter = [];
    end

    qv_max_full = zeros(mpc.nh, mpc.N);
    if ~isempty(qv_max)
        qv_max_full = fill_vec(qv_max_full, qv_max, 1);
    end
    h_cnstr.qv_max = qv_max_full(:,1:mpc.N-1);
    if ~isempty(h_cnstr.use_k0)
        h_cnstr.qv_max_0 = qv_max_full(h_cnstr.rows_k0,1);
    else
        h_cnstr.qv_max_0 = [];
    end
    if ~isempty(h_cnstr.use_ter)
        h_cnstr.qv_max_ter = qv_max_full(h_cnstr.rows_ter,mpc.N);
    else
        h_cnstr.qv_max_ter = [];
    end

else
    h_cnstr.max_limit = [];
    h_cnstr.max = [];
    h_cnstr.max_0 = [];
    h_cnstr.max_ter = [];
    h_cnstr.qv_max = [];
    h_cnstr.qv_max_0 = [];
    h_cnstr.qv_max_ter = [];
end

mpc.h_cnstr = h_cnstr;

end
