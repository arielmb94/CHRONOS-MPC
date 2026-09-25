% INIT_MPC_CUSTOM_COST Add a quadratic or linear cost on a user-defined signal.
%
%   mpc = INIT_MPC_CUSTOM_COST(mpc, Cz, Dz, Dsuz, Ddz, Qz, qz) defines
%
%       z_k = Cz_k*s_k + Dz_k*u_k + Dsuz_k*su_k + Ddz_k*dz_k
%
%   and adds quadratic and/or linear penalties on z_k:
%
%       J_custom += 0.5*z_k'*Qz_k*z_k + qz_k'*z_k
%
%   Use [] for a coefficient or penalty that is not needed. Here, su_k is
%   the control action preceding u_k, and dz_k is a dedicated fixed known
%   input supplied through the dz_in argument of MPC_SOLVE.
%
%   Coefficient and quadratic-weight matrices may be constant or contain L
%   horizon stages in their third dimension. The linear weight qz may be an
%   nz-by-1 vector or an nz-by-L matrix. L is the number of supplied horizon
%   stages; if L < N, the last supplied stage is reused for the remaining
%   stages.
%
%   Call this function after INIT_MPC_DYNAMICS and before
%   BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure.
%     Cz      - State coefficient, size nz-by-nx or nz-by-nx-by-L.
%     Dz      - Control coefficient, size nz-by-nu or nz-by-nu-by-L.
%     Dsuz    - Previous-control coefficient, size nz-by-nu or
%               nz-by-nu-by-L.
%     Ddz     - Fixed-known-input coefficient, size nz-by-ndz or
%               nz-by-ndz-by-L.
%     Qz      - Optional quadratic weight, size nz-by-nz or nz-by-nz-by-L.
%     qz      - Optional linear weight, size nz-by-1 or nz-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - penalize deviation from a runtime control reference:
%
%       % z = u - u_ref
%       Cz   = zeros(nu, nx);
%       Dz   = eye(nu);
%       Dsuz = [];
%       Ddz  = -eye(nu);
%       Qz   = eye(nu);
%       mpc = init_mpc_Custom_cost(mpc, Cz, Dz, Dsuz, Ddz, Qz, []);
%
%   Pass u_ref as the dz_in argument of MPC_SOLVE at runtime.
function mpc = init_mpc_Custom_cost(mpc,Cz,Dz,Dsuz,Ddz,Qz,qz)
arguments
    mpc
    Cz = []
    Dz = []
    Dsuz = []
    Ddz = []
    Qz = [];
    qz = [];
end

% Empty and all-zero optional terms are omitted before dimension inference or
% validation, independently of their supplied dimensions.
if ~isempty(Cz) && ~any(Cz(:)), Cz = []; end
if ~isempty(Dz) && ~any(Dz(:)), Dz = []; end
if ~isempty(Dsuz) && ~any(Dsuz(:)), Dsuz = []; end
if ~isempty(Ddz) && ~any(Ddz(:)), Ddz = []; end
if ~isempty(Qz) && ~any(Qz(:)), Qz = []; end
if ~isempty(qz) && ~any(qz(:)), qz = []; end

% A custom cost needs both a decision-dependent signal and an active weight.
if (isempty(Cz) && isempty(Dz) && isempty(Dsuz)) || ...
        (isempty(Qz) && isempty(qz))
    return;
end

%number of custom costs
if ~isempty(Cz)
    mpc.nz = size(Cz,1);  
elseif ~isempty(Dz)
    mpc.nz = size(Dz,1);
else
    mpc.nz = size(Dsuz,1);
end

validate_matrix(Cz, mpc.nz, mpc.nx, 'Cz');
validate_matrix(Dz, mpc.nz, mpc.nu, 'Dz');
validate_matrix(Dsuz, mpc.nz, mpc.nu, 'Dsuz');
if ~isempty(Ddz)
    mpc.ndz = size(Ddz,2);
    validate_matrix(Ddz, mpc.nz, mpc.ndz, 'Ddz');
end

validate_matrix(Qz, mpc.nz, mpc.nz, 'Qz', true);
validate_column_vector(qz, mpc.nz, 'qz', true);

if isscalar(Qz)
    Qz = Qz * eye(mpc.nz);
elseif size(Qz,1) == 1 && size(Qz,2) == 1
    Qz_staged = Qz;
    Qz = zeros(mpc.nz, mpc.nz, size(Qz_staged,3));
    Qz_eye = eye(mpc.nz);
    for k = 1:size(Qz_staged,3)
        Qz(:,:,k) = Qz_staged(1,1,k) * Qz_eye;
    end
end

if ~isempty(qz) && size(qz,1) == 1
    qz = ones(mpc.nz, 1) * qz;
end

if ~isempty(Cz)
    mpc.z_use_s = 1;

    len_Cz = size(Cz,3);
    if len_Cz < mpc.N
        mpc.Cz = zeros(mpc.nz,mpc.nx,mpc.N-1);
        mpc.Cz = fill_mat(mpc.Cz, Cz, 1);
        Cz_ter = mpc.Cz(:,:,mpc.N-1);
    else
        mpc.Cz = Cz(:,:,1:mpc.N-1);
        Cz_ter = Cz(:,:,mpc.N);
    end
end

if ~isempty(Dz)
    mpc.z_use_u = 1;

    len_Dz = size(Dz,3);
    if len_Dz < mpc.N
        mpc.Dz = zeros(mpc.nz,mpc.nu,mpc.N-1);
        mpc.Dz = fill_mat(mpc.Dz, Dz, 1);
        Dz_ter = mpc.Dz(:,:,mpc.N-1);
    else
        mpc.Dz = Dz(:,:,1:mpc.N-1);
        Dz_ter = Dz(:,:,mpc.N);
    end
end

if ~isempty(Dsuz)
    mpc.z_use_su = 1;

    mpc.has_du = 1;

    len_Dsuz = size(Dsuz,3);
    if len_Dsuz < mpc.N
        mpc.Dsuz = zeros(mpc.nz,mpc.nu,mpc.N-1);
        mpc.Dsuz = fill_mat(mpc.Dsuz, Dsuz, 1);
        Dsuz_ter = mpc.Dsuz(:,:,mpc.N-1);
    else
        mpc.Dsuz = Dsuz(:,:,1:mpc.N-1);
        Dsuz_ter = Dsuz(:,:,mpc.N);
    end
end

if ~isempty(Ddz)
    mpc.z_use_d = 1;

    len_Ddz = size(Ddz,3);
    if len_Ddz < mpc.N
        mpc.Ddz = zeros(mpc.nz,mpc.ndz,mpc.N-1);
        mpc.Ddz = fill_mat(mpc.Ddz, Ddz, 1);
        Ddz_ter = mpc.Ddz(:,:,mpc.N-1);
    else
        mpc.Ddz = Ddz(:,:,1:mpc.N-1);
        Ddz_ter = Ddz(:,:,mpc.N);
    end
end


mpc.z_use_k0 = [];
mpc.z_use_ter = [];
mpc.Cz_0 = [];
mpc.Dz_0 = [];
mpc.Dsuz_0 = [];
mpc.Ddz_0 = [];
mpc.Cz_ter = [];
mpc.z_rows_k0 = [];
mpc.z_rows_ter = [];
% at k = 0, only rows with Dz!=0 (with dependence on control action u) are
% considered
z_row_0 = find(~all(Dz(:,:,1)==0,2));
mpc.nz_0 = length(z_row_0);

if mpc.nz_0
    mpc.z_rows_k0 = z_row_0;
    mpc.z_use_k0 = 1;

    if ~isempty(mpc.z_use_s), mpc.Cz_0 = Cz(mpc.z_rows_k0,:,1); end
    if ~isempty(mpc.z_use_u), mpc.Dz_0 = Dz(mpc.z_rows_k0,:,1); end
    if ~isempty(mpc.z_use_su), mpc.Dsuz_0 = Dsuz(mpc.z_rows_k0,:,1); end
    if ~isempty(mpc.z_use_d), mpc.Ddz_0 = Ddz(mpc.z_rows_k0,:,1); end
end

% at k = N, only rows strictly dependent on s are considered
if ~isempty(mpc.z_use_s)
    strict_s_rows = any(Cz_ter~=0,2);
    if ~isempty(mpc.z_use_u), strict_s_rows = strict_s_rows & all(Dz_ter==0,2); end
    if ~isempty(mpc.z_use_su), strict_s_rows = strict_s_rows & all(Dsuz_ter==0,2); end
    if ~isempty(mpc.z_use_d), strict_s_rows = strict_s_rows & all(Ddz_ter==0,2); end

    z_row_ter = find(strict_s_rows);
    mpc.nz_ter = length(z_row_ter);
else
    mpc.nz_ter = 0;
end
if mpc.nz_ter
    mpc.z_rows_ter = z_row_ter;
    mpc.z_use_ter = 1;

    mpc.Cz_ter = Cz_ter(mpc.z_rows_ter,:); 
end

% init z vector
if ~isempty(mpc.z_use_k0), mpc.z_0 = zeros(mpc.nz_0,1); else, mpc.z_0=[]; end
mpc.z = zeros(mpc.nz,mpc.N-1);
if ~isempty(mpc.z_use_ter), mpc.z_ter = zeros(mpc.nz_ter,1); else, mpc.z_ter=[]; end
% init d vector
if mpc.ndz
    mpc.dz = zeros(mpc.ndz,mpc.N);
end

% Quadratic cost matrix
if any(Qz(:))
    mpc.quad_custom_cost = 1;

    mpc.Qz = zeros(mpc.nz,mpc.nz,mpc.N-1);
    len_Qz = size(Qz,3);
    if len_Qz < mpc.N
        mpc.Qz = fill_mat(mpc.Qz, Qz, 1);
        Qz_ter = mpc.Qz(:,:,mpc.N-1);
    else
        mpc.Qz = Qz(:,:,1:mpc.N-1);
        Qz_ter = Qz(:,:,mpc.N);
    end

    if ~isempty(mpc.z_use_k0)
        mpc.Qz_0 = Qz(mpc.z_rows_k0,mpc.z_rows_k0,1);
    else
        mpc.Qz_0 = [];
    end
    if ~isempty(mpc.z_use_ter)
        mpc.Qz_ter = Qz_ter(mpc.z_rows_ter,mpc.z_rows_ter);
    else
        mpc.Qz_ter = [];
    end
else
    mpc.Qz_0 = [];
    mpc.Qz_ter = [];
end

% Linear cost matrix
if any(qz(:))
    mpc.lin_custom_cost = 1;

    mpc.qz = zeros(mpc.nz,mpc.N-1);
    len_qz = size(qz,2);
    if len_qz < mpc.N
        mpc.qz = fill_vec(mpc.qz, qz, 1);
        qz_ter = mpc.qz(:,mpc.N-1);
    else
        mpc.qz = qz(:,1:mpc.N-1);
        qz_ter = qz(:,mpc.N);
    end

    if ~isempty(mpc.z_use_k0)
        mpc.qz_0 = qz(mpc.z_rows_k0,1);
    else
        mpc.qz_0 = [];
    end
    if ~isempty(mpc.z_use_ter)
        mpc.qz_ter = qz_ter(mpc.z_rows_ter);
    else
        mpc.qz_ter = [];
    end
else
    mpc.qz_0 = [];
    mpc.qz_ter = [];
end
 
end
