% INIT_MPC_DYNAMICS Define the prediction model.
%
%   mpc = INIT_MPC_DYNAMICS(mpc, A, B) defines
%
%       s_(k+1) = A_k*s_k + B_k*u_k.
%
%   mpc = INIT_MPC_DYNAMICS(mpc, A, B, Bd) adds the fixed known input d_k
%
%       s_(k+1) = A_k*s_k + B_k*u_k + Bd_k*d_k.
% 
%   Pass d_k as the d_in argument when building and solving the
%   MPC. Use [] for Bd when the model has no such input.
%
%   A, B, and Bd may be constant matrices or three-dimensional arrays. For
%   an array, L is the number of supplied horizon stages. If L < N, the last
%   stage is reused through the rest of the horizon; if L >= N, the first N
%   stages are used.
%
%   INIT_MPC_DYNAMICS automatically sets the tracking output to the full state,
%   y_k = I*s_k = s_k. If this is the desired output, no separate output
%   initialization is needed. Call INIT_MPC_OUTPUT after this function only
%   when a different tracking or constrained output is required.
%
%   Call this function after INIT_MPC and before adding costs or
%   constraints.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure created by INIT_MPC.
%     A       - State matrix, size nx-by-nx or nx-by-nx-by-L.
%     B       - Control matrix, size nx-by-nu or nx-by-nu-by-L.
%     Bd      - Optional fixed-known-input matrix, size nx-by-nd or
%               nx-by-nd-by-L. Default: [].
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - define a time-invariant model without a known input:
%
%       mpc = init_mpc(20);
%       mpc = init_mpc_dynamics(mpc, A, B, []);
function mpc = init_mpc_dynamics(mpc,A,B,Bd)
arguments
    mpc
    A = []
    B = []
    Bd = []
end

% An all-zero optional disturbance matrix means that the term is omitted.
if ~isempty(Bd) && ~any(Bd(:))
    Bd = [];
end

% Infer and validate all model dimensions before allocating MPC fields.
nx = size(A,1);
validate_matrix(A, nx, nx, 'A');

nu = size(B,2);
validate_matrix(B, nx, nu, 'B');

if ~isempty(Bd)
    nd = size(Bd,2);
    validate_matrix(Bd, nx, nd, 'Bd');
end

mpc.Bd = Bd;

%number of states
mpc.nx = nx;
mpc.A = zeros(mpc.nx,mpc.nx,mpc.N);
mpc.A = fill_mat(mpc.A, A, 1);
%number of control inputs
mpc.nu = nu;
mpc.B = zeros(mpc.nx,mpc.nu,mpc.N);
mpc.B = fill_mat(mpc.B, B, 1);
%number of disturbance inputs
if ~isempty(Bd)
    mpc.nd = nd;
    mpc.dyn_use_d = 1;
    mpc.Bd = zeros(mpc.nx,mpc.nd,mpc.N);
    mpc.Bd = fill_mat(mpc.Bd, Bd, 1);
end

mpc.s = zeros(mpc.nx,mpc.N);
mpc.s_ter = zeros(mpc.nx,1);
mpc.su = zeros(mpc.nu,mpc.N);
mpc.u = zeros(mpc.nu,mpc.N);
mpc.du = zeros(mpc.nu,mpc.N);

if ~isempty(mpc.dyn_use_d)
    mpc.d = zeros(mpc.nd,mpc.N);
end

% Assume C = I*x
mpc.D = [];
mpc.Dd = [];
mpc.C_0 = [];
mpc.D_0 = [];
mpc.Dd_0 = [];

C = eye(mpc.nx);
mpc.C = zeros(mpc.nx,mpc.nx,mpc.N-1);
mpc.C = fill_mat(mpc.C, C, 1);
mpc.C_ter = eye(mpc.nx);

mpc.ny = mpc.nx;
mpc.ny_0 = 0;
mpc.ny_ter = mpc.nx;

mpc.y_use_k0 = [];
mpc.y_rows_k0 = [];
mpc.y_use_ter = 1;
mpc.y_rows_ter = 1:mpc.nx; 

mpc.y_use_s = 1;
mpc.y_use_u = [];
mpc.y_use_d = [];

% init y, reference and error vectors
mpc.r_0 = zeros(mpc.ny_0,1);
mpc.y_0 = zeros(mpc.ny_0,1);
mpc.err_0 = zeros(mpc.ny_0,1);
mpc.r = zeros(mpc.ny,mpc.N-1);
mpc.y = zeros(mpc.ny,mpc.N-1);
mpc.err = zeros(mpc.ny,mpc.N-1);
mpc.r_ter = zeros(mpc.ny_ter,1);
mpc.y_ter = zeros(mpc.ny_ter,1);
mpc.err_ter = zeros(mpc.ny_ter,1);

end
