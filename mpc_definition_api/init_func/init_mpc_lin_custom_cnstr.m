% INIT_MPC_LIN_CUSTOM_CNSTR Defines constraints and soft-penalties on
% custom user defined signals.
%
%   mpc = INIT_MPC_LIN_CUSTOM_CNSTR(mpc, Ch, Dh, Ddh, h_min, h_max) defines a 
%   custom auxiliary signal 'h' and sets strict (hard) lower and upper bounds 
%   on it. The custom signal is calculated as:
%
%       h = Ch * x + Dh * u + Ddh * dh
%
%   The solver will strictly enforce h_min <= h <= h_max. This is suitable for 
%   hard physical limits, but massive disturbances could cause the solver to crash 
%   (go infeasible) if the limits are impossible to satisfy.
%
%   mpc = INIT_MPC_LIN_CUSTOM_CNSTR(mpc, ..., h_min_slack_active, h_max_slack_active, qv_min, qv_max) 
%   allows you to define these custom bounds as "soft" constraints. Soft constraints 
%   can be safely violated during severe disturbances to prevent solver crashes, 
%   while applying a customizable linear penalty to drive the signal back within 
%   limits as quickly as actuator power allows. It is highly advisible to
%   enable soft constraint on the custom signals.
%
%   INPUTS:
%       mpc                - CHRONOS MPC structure.
%       Ch                 - [nh x nx] Matrix mapping states to the custom signal.
%       Dh                 - [nh x nu] Matrix mapping inputs to the custom signal.
%       Ddh                - [nh x ndh] Matrix mapping measured disturbances to the custom signal.
%       h_min              - [nh x 1] Array of lower limits (use [] if none).
%       h_max              - [nh x 1] Array of upper limits (use [] if none).
%       qv_min             - (Optional) [nh x 1] or scalar. Penalty weight for violating 
%                            the h_min soft limits. Higher values mean stricter enforcement.
%       qv_max             - (Optional) [nh x 1] or scalar. Penalty weight for violating 
%                            the h_max soft limits. Higher values mean stricter enforcement.
%
%   OUTPUTS:
%       mpc                - Updated MPC structure. All necessary background math 
%                            (matrices, gradients, Hessians, and slack variables) 
%                            are automatically assembled and added to the object.
%
%   EXAMPLE USE CASE: Tracking an input reference (u_star)
%       We want to limit the variation of the control action with respect to a 
%       dynamic target value, meaning we want to constrain: h = u - u_star.
%       To achieve this, we set up our custom signal matrices as:
%           - Ch  = zeros(nu, nx)
%           - Dh  = eye(nu)
%           - Ddh = -eye(nu)
%       This creates the equation: h = 0*x + 1*u - 1*dh.
%       During runtime, the user passes 'u_star' into the 'dh' disturbance 
%       vector, and the solver handles the rest!
%
%   USAGE TIPS:
%       - If qv_min or qv_max are not passed, the soft constraint penalty 
%         weight will default to the value stored in mpc.qv
%       - Passing a scalar to the slack or qv inputs will automatically apply 
%         that setting across all constrained outputs.
function mpc = init_mpc_lin_custom_cnstr(mpc,Ch,Dh,Dsuh,Ddh,...
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

% general constraints boolean
mpc.has_h_cnstr = 1;

%number of general inequalities
if ~isempty(Ch) && any(Ch(:))
    mpc.nh = size(Ch,1);  
elseif ~isempty(Dh) && any(Dh(:))
    mpc.nh = size(Dh,1);
elseif ~isempty(Dsuh) && any(Dsuh(:))
    mpc.nh = size(Dsuh,1);
end

% INPUT DIMENSION VALIDATION 
validate_column_vector(h_min, mpc.nh, 'h_min');
validate_column_vector(h_max, mpc.nh, 'h_max');
validate_column_vector(qv_min, mpc.nh, 'qv_min');
validate_column_vector(qv_max, mpc.nh, 'qv_max');

h_cnstr.use_s = 0;
h_cnstr.use_u = 0;
h_cnstr.use_su = 0;
h_cnstr.use_d = 0;

if ~isempty(Ch) && any(Ch(:))
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

if ~isempty(Dh) && any(Dh(:))
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

if ~isempty(Dsuh) && any(Dsuh(:))
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

if ~isempty(Ddh) && any(Ddh(:))
    h_cnstr.use_d = 1;

    mpc.ndh = size(Ddh,2);

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

h_cnstr.use_k0 = 0;
h_cnstr.use_ter = 0;
% at k = 0, only rows with Dh!=0 (with dependence on control action u) are
% considered
h_row_0 = find(~all(Dh(:,:,1)==0,2));
mpc.nh_0 = length(h_row_0);

if mpc.nh_0
    h_cnstr.rows_k0 = h_row_0;
    h_cnstr.use_k0 = 1;

    if h_cnstr.use_s, mpc.Ch_0 = Ch(h_cnstr.rows_k0,:,1); end
    if h_cnstr.use_u, mpc.Dh_0 = Dh(h_cnstr.rows_k0,:,1); end
    if h_cnstr.use_su, mpc.Dsuh_0 = Dsuh(h_cnstr.rows_k0,:,1); end
    if h_cnstr.use_d, mpc.Ddh_0 = Ddh(h_cnstr.rows_k0,:,1); end
end

% at k = N, only rows strictly dependent on s are considered
if  h_cnstr.use_s
    strict_s_rows = any(Ch_ter~=0,2);
    if h_cnstr.use_u, strict_s_rows = strict_s_rows & all(Dh_ter==0,2); end
    if h_cnstr.use_su, strict_s_rows = strict_s_rows & all(Dsuh_ter==0,2); end
    if h_cnstr.use_d, strict_s_rows = strict_s_rows & all(Ddh_ter==0,2); end

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
if h_cnstr.use_k0, mpc.h_0 = zeros(mpc.nh_0,1); else, mpc.h_0=[]; end
mpc.h = zeros(mpc.nh,mpc.N-1);
if h_cnstr.use_ter, mpc.h_ter = zeros(mpc.nh_ter,1); else, mpc.h_ter=[]; end

% init dh disturbance vector
if h_cnstr.use_d
    mpc.dh = zeros(mpc.ndh,mpc.N);
end

% Expand scalars to full vectors if needed
if isscalar(h_min), h_min = h_min * ones(mpc.nh, 1); end
if isscalar(h_max), h_max = h_max * ones(mpc.nh, 1); end
if isscalar(qv_min), qv_min = qv_min * ones(mpc.nh, 1); end
if isscalar(qv_max), qv_max = qv_max * ones(mpc.nh, 1); end

if ~isempty(h_min)

    h_cnstr.min_limit = 1;

    if h_cnstr.use_k0
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.nh_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.nh_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.nh;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.nh;
    if h_cnstr.use_ter
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.nh_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.nh_ter;
    end

    h_cnstr.g_min_index_k = [];
    h_cnstr.v_min_index_k = [];

    h_min_full = zeros(mpc.nh, mpc.N);
    h_min_full = fill_vec(h_min_full, h_min, 1);
    h_cnstr.min = h_min_full(:,1:mpc.N-1);
    if h_cnstr.use_k0, h_cnstr.min_0 = h_min_full(h_cnstr.rows_k0,1); end
    if h_cnstr.use_ter, h_cnstr.min_ter = h_min_full(h_cnstr.rows_ter,mpc.N); end

    qv_min_full = zeros(mpc.nh, mpc.N);
    if ~isempty(qv_min)
        qv_min_full = fill_vec(qv_min_full, qv_min, 1);
    end
    h_cnstr.qv_min = qv_min_full(:,1:mpc.N-1);
    if h_cnstr.use_k0, h_cnstr.qv_min_0 = qv_min_full(h_cnstr.rows_k0,1); end
    if h_cnstr.use_ter, h_cnstr.qv_min_ter = qv_min_full(h_cnstr.rows_ter,mpc.N); end
   
else
    h_cnstr.min_limit = 0;
end

if ~isempty(h_max)

    h_cnstr.max_limit = 1;

    if h_cnstr.use_k0
        mpc.ng_k(1) = mpc.ng_k(1) + mpc.nh_0;
        mpc.nv_k(1) = mpc.nv_k(1) + mpc.nh_0;
    end
    mpc.ng_k(2) = mpc.ng_k(2) + mpc.nh;
    mpc.nv_k(2) = mpc.nv_k(2) + mpc.nh;
    if h_cnstr.use_ter
        mpc.ng_k(3) = mpc.ng_k(3) + mpc.nh_ter;
        mpc.nv_k(3) = mpc.nv_k(3) + mpc.nh_ter;
    end

    h_cnstr.g_max_index_k = [];
    h_cnstr.v_max_index_k = [];

    h_max_full = zeros(mpc.nh, mpc.N);
    h_max_full = fill_vec(h_max_full, h_max, 1);
    h_cnstr.max = h_max_full(:,1:mpc.N-1);
    if h_cnstr.use_k0, h_cnstr.max_0 = h_max_full(h_cnstr.rows_k0,1); end
    if h_cnstr.use_ter, h_cnstr.max_ter = h_max_full(h_cnstr.rows_ter,mpc.N); end

    qv_max_full = zeros(mpc.nh, mpc.N);
    if ~isempty(qv_max)
        qv_max_full = fill_vec(qv_max_full, qv_max, 1);
    end
    h_cnstr.qv_max = qv_max_full(:,1:mpc.N-1);
    if h_cnstr.use_k0, h_cnstr.qv_max_0 = qv_max_full(h_cnstr.rows_k0,1); end
    if h_cnstr.use_ter, h_cnstr.qv_max_ter = qv_max_full(h_cnstr.rows_ter,mpc.N); end

else
    h_cnstr.max_limit = 0;
end

mpc.h_cnstr = h_cnstr;

end