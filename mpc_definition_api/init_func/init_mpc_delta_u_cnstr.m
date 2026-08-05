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

mpc.has_du_cnstr = 1;
mpc.has_du = 1;

% INPUT DIMENSION VALIDATION 
validate_column_vector(du_min, mpc.nu, 'du_min');
validate_column_vector(du_max, mpc.nu, 'du_max');

% Expand scalars to full local vectors if needed
if isscalar(du_min), du_min = du_min * ones(mpc.nu, 1); end
if isscalar(du_max), du_max = du_max * ones(mpc.nu, 1); end

du_cnstr.use_k0 = 0;
du_cnstr.use_ter = 0;
du_cnstr.rows_k0 = [];
du_cnstr.rows_ter = [];

if ~isempty(du_min)

    du_cnstr.min_limit = 1;
    
    mpc.ng_k(1:2) = mpc.ng_k(1:2) + mpc.nu;

    du_cnstr.g_min_index_k = [];

    du_cnstr.min = zeros(mpc.nu,mpc.N);
    du_cnstr.min = fill_vec(du_cnstr.min, du_min, 1);
else
    du_cnstr.min_limit = 0;
end

if ~isempty(du_max)

    du_cnstr.max_limit = 1;

    mpc.ng_k(1:2) = mpc.ng_k(1:2) + mpc.nu;

    du_cnstr.g_max_index_k = [];

    du_cnstr.max = zeros(mpc.nu,mpc.N);
    du_cnstr.max = fill_vec(du_cnstr.max, du_max, 1);    
else
    du_cnstr.max_limit = 0;
end

mpc.du_cnstr = du_cnstr;

end