%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Ddz,Qz,qz)
%
% Adds quadratic and linear penalties on a custom user defined signal z:
%
%   J += z'*Qz*z + qz*z
%
% The signal z is defined as:
%
%   z = Cz * x + Dz * u + Ddz * dz
%
% The user must define the signal z by selecting appropiate values for the
% matrices Cz, Dz, Ddz.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - Cz: nz x nx matrix, states output matrix
%   - Dz: nz x nu matrix, input feedtrhough matrix
%   - Ddz: nz x ndi matrix, disturbance feedtrhough matrix
%   - Qz (optional): nz x nz square matrix, weights for the quadratic
%   penalty term on the user defined signal z
%   - qz (optional): nz column vector, weights for the linear penalty term
%   on the user defined signal z. Use the linear penalty term only in the
%   case that z takes only positive values
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
% Example:
% We want to penalize the variation of the MPC control action with respect 
% a given value u_star, e.g. we want to minimize z = u - u_star.
% To achieve this we define z by selecting:
%   - Cz = [0 0 ... 0]
%   - Dz = [1]
%   - Ddz = [-1]
% which corresponds to z = [0 0 ... 0] * x + [1] * u + [-1] * di
% The "disturbance" term on z corresponds to u_star, to be introduced on 
% the appropiate field on mpc_solve() during runtime MPC execution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = init_mpc_Lin_Custom_cost(mpc,Cz,Dz,Dsuz,Ddz,Qz,qz)
arguments
    mpc
    Cz = []
    Dz = []
    Dsuz = []
    Ddz = []
    Qz = [];
    qz = [];
end

%number of custom costs
if ~isempty(Cz) && any(Cz(:))
    mpc.nz = size(Cz,1);  
elseif ~isempty(Dz) && any(Dz(:))
    mpc.nz = size(Dz,1);
elseif ~isempty(Dsuz) && any(Dsuz(:))
    mpc.nz = size(Dsuz,1);
end

if ~isempty(Cz) && any(Cz(:))
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

if ~isempty(Dz) && any(Dz(:))
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

if ~isempty(Dsuz) && any(Dsuz(:))
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

if ~isempty(Ddz) && any(Ddz(:))
    mpc.z_use_d = 1;

    mpc.ndz = size(Ddz,2);

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


mpc.z_use_k0 = 0;
mpc.z_use_ter = 0;
% at k = 0, only rows with Dz!=0 (with dependence on control action u) are
% considered
z_row_0 = find(~all(Dz(:,:,1)==0,2));
mpc.nz_0 = length(z_row_0);

if mpc.nz_0
    mpc.z_rows_k0 = z_row_0;
    mpc.z_use_k0 = 1;

    if mpc.z_use_s, mpc.Cz_0 = Cz(mpc.z_rows_k0,:,1); end
    if mpc.z_use_u, mpc.Dz_0 = Dz(mpc.z_rows_k0,:,1); end
    if mpc.z_use_su, mpc.Dsuz_0 = Dsuz(mpc.z_rows_k0,:,1); end
    if mpc.z_use_d, mpc.Ddz_0 = Ddz(mpc.z_rows_k0,:,1); end
end

% at k = N, only rows strictly dependent on s are considered
if  mpc.z_use_s
    strict_s_rows = any(Cz_ter~=0,2);
    if mpc.z_use_u, strict_s_rows = strict_s_rows & all(Dz_ter==0,2); end
    if mpc.z_use_su, strict_s_rows = strict_s_rows & all(Dsuz_ter==0,2); end
    if mpc.z_use_d, strict_s_rows = strict_s_rows & all(Ddz_ter==0,2); end

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
if mpc.z_use_k0, mpc.z_0 = zeros(mpc.nz_0,1); else, mpc.z_0=[]; end
mpc.z = zeros(mpc.nz,mpc.N-1);
if mpc.z_use_ter, mpc.z_ter = zeros(mpc.nz_ter,1); else, mpc.z_ter=[]; end
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

    if mpc.z_use_k0, mpc.Qz_0 = Qz(mpc.z_rows_k0,mpc.z_rows_k0,1); end
    if mpc.z_use_ter, mpc.Qz_ter = Qz_ter(mpc.z_rows_ter,mpc.z_rows_ter); end
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

    if mpc.z_use_k0, mpc.qz_0 = qz(mpc.z_rows_k0,1); end
    if mpc.z_use_ter, mpc.qz_ter = qz_ter(mpc.z_rows_ter); end
end
 
end