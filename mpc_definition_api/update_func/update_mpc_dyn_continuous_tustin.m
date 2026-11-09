function mpc = update_mpc_dyn_continuous_tustin(mpc, Ac, Bc, Bdc, Ts, M)
% UPDATE_MPC_DYN_CONTINUOUS_TUSTIN Update dynamics with the Tustin method.
%
%   mpc = UPDATE_MPC_DYN_CONTINUOUS_TUSTIN(mpc, Ac, Bc, Bdc, Ts) updates
%   the full prediction model by discretizing
%
%       x_dot = Ac*x + Bc*u + Bdc*d.
%
%   mpc = UPDATE_MPC_DYN_CONTINUOUS_TUSTIN(mpc, Ac, Bc, Bdc, Ts, M) uses
%   the optional mass formulation
%
%       M*x_dot = Ac*x + Bc*u + Bdc*d.
%
%   Pass the continuous coefficients directly; discretization incorporates
%   M. Omitted, empty, or scalar-zero M means I. The Tustin model is
%
%       (M - (Ts/2)*Ac)*x_(k+1) = (M + (Ts/2)*Ac)*x_k
%                                + Ts*Bc*u_k + Ts*Bdc*d_k.
%
%   Pass full nonempty Ac and Bc on every call. If the known-input matrix Bd
%   was defined at initialization, also pass full nonempty Bdc. Otherwise
%   Bdc is ignored; use [].
%   Correctly sized zero coefficients update numerical data, including
%   zeroing an existing Bd. The built dimensions and enabled terms stay fixed.
%
%   Ac, Bc, Bdc, and M may independently be constant matrices or contain L
%   pages. Each short input repeats its last page through horizon N; pages
%   beyond N are ignored. Unused Bdc pages do not affect staging.
%
%   Ts must be a positive scalar. Nondefault M must be invertible.
%   Inputs must have the built dimensions; no runtime dimension or numerical
%   checks are performed.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure with initialized dynamics.
%     Ac      - State coefficient, size nx-by-nx or nx-by-nx-by-L.
%     Bc      - Control coefficient, size nx-by-nu or nx-by-nu-by-L.
%     Bdc     - Required if known-input matrix Bd was defined at
%               initialization, size nx-by-nd or nx-by-nd-by-L; ignored
%               otherwise.
%     Ts      - Positive scalar sampling interval.
%     M       - Optional mass matrix, size nx-by-nx or nx-by-nx-by-L.
%               Default: I for omitted, empty, or scalar-zero M.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure; retain this returned value.
%
%   Example - update a model built without a known input:
%
%       mpc = update_mpc_dyn_continuous_tustin(mpc, Ac, Bc, [], Ts);
%
%   Example - update staged state and mass matrices with a known input:
%
%       Ac = cat(3, Ac_1, Ac_2);
%       M = cat(3, M_1, M_2, M_3);
%       mpc = update_mpc_dyn_continuous_tustin(mpc, Ac, Bc, Bdc, Ts, M);

nx = mpc.nx;
N = mpc.N;
has_dist = ~isempty(mpc.dyn_use_d);
has_M = (nargin >= 6) && ~isempty(M) && ~(isscalar(M) && M == 0);

len_Ac = size(Ac,3);
len_Bc = size(Bc,3);
len_Bdc = 0;
if has_dist
    len_Bdc = size(Bdc,3);
end
len_M = 0;
if has_M
    len_M = size(M,3);
end
num_stages = min(N, max([len_Ac, len_Bc, len_Bdc, len_M]));

% Inherit fixed stage shapes from the built controller.
Ac_k = mpc.A(:,:,1);
Bc_k = mpc.B(:,:,1);
if has_dist
    Bdc_k = mpc.Bd(:,:,1);
end
M_k = mpc.A(:,:,1);
W = mpc.A(:,:,1);
rhs_A = mpc.A(:,:,1);
if ~has_M
    M_k(:,:) = 0;
    for j = 1:nx
        M_k(j,j) = 1;
    end
end

for k = 1:num_stages
    Ac_k(:,:) = Ac(:,:,min(k,len_Ac));
    Bc_k(:,:) = Bc(:,:,min(k,len_Bc));
    if has_dist
        Bdc_k(:,:) = Bdc(:,:,min(k,len_Bdc));
    end
    if has_M
        M_k(:,:) = M(:,:,min(k,len_M));
    end

    W(:,:) = M_k - (Ts/2) * Ac_k;
    rhs_A(:,:) = M_k + (Ts/2) * Ac_k;
    Ac_k(:,:) = solveLinearSystemLU(W, rhs_A);
    Bc_k(:,:) = solveLinearSystemLU(W, Ts * Bc_k);
    if has_dist
        Bdc_k(:,:) = solveLinearSystemLU(W, Ts * Bdc_k);
    end

    mpc.A(:,:,k) = Ac_k;
    mpc.B(:,:,k) = Bc_k;
    if has_dist
        mpc.Bd(:,:,k) = Bdc_k;
    end
end

% Reuse the final discrete stage without repeating its conversion.
for k = num_stages+1:N
    mpc.A(:,:,k) = mpc.A(:,:,num_stages);
    mpc.B(:,:,k) = mpc.B(:,:,num_stages);
    if has_dist
        mpc.Bd(:,:,k) = mpc.Bd(:,:,num_stages);
    end
end
end
