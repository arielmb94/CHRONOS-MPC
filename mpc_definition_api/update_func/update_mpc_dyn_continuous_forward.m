function mpc = update_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts, M)
% UPDATE_MPC_DYN_CONTINUOUS_FORWARD Update dynamics with forward Euler.
%
%   mpc = UPDATE_MPC_DYN_CONTINUOUS_FORWARD(mpc, Ac, Bc, Bdc, Ts) updates
%   the full prediction model by discretizing
%
%       x_dot = Ac*x + Bc*u + Bdc*d.
%
%   mpc = UPDATE_MPC_DYN_CONTINUOUS_FORWARD(mpc, Ac, Bc, Bdc, Ts, M) uses
%   the optional mass formulation
%
%       M*x_dot = Ac*x + Bc*u + Bdc*d.
%
%   Pass the continuous coefficients directly; discretization incorporates
%   M. Omitted, empty, or scalar-zero M means I. Forward Euler gives
%
%       A = I + Ts*M^(-1)*Ac, B = Ts*M^(-1)*Bc,
%       Bd = Ts*M^(-1)*Bdc.
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
%       mpc = update_mpc_dyn_continuous_forward(mpc, Ac, Bc, [], Ts);
%
%   Example - update staged state and mass matrices with a known input:
%
%       Ac = cat(3, Ac_1, Ac_2);
%       M = cat(3, M_1, M_2, M_3);
%       mpc = update_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts, M);

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
I = mpc.A(:,:,1);
I(:,:) = 0;
for j = 1:nx
    I(j,j) = 1;
end
Ac_k = mpc.A(:,:,1);
Bc_k = mpc.B(:,:,1);
if has_dist
    Bdc_k = mpc.Bd(:,:,1);
end
M_k = mpc.A(:,:,1);
L = mpc.A(:,:,1);
U = mpc.A(:,:,1);

for k = 1:num_stages
    first_stage = (k == 1);
    factor_changed = has_M && (first_stage || k <= len_M);
    if factor_changed
        M_k(:,:) = M(:,:,min(k,len_M));
        [L(:,:), U(:,:)] = factorizeLinearSystemLU(M_k);
    end

    if first_stage || k <= len_Ac || factor_changed
        Ac_k(:,:) = Ac(:,:,min(k,len_Ac));
        if has_M
            Ac_k(:,:) = I + Ts * solveFactoredLinearSystemLU(L, U, Ac_k);
        else
            Ac_k(:,:) = I + Ts * Ac_k;
        end
    end
    if first_stage || k <= len_Bc || factor_changed
        Bc_k(:,:) = Bc(:,:,min(k,len_Bc));
        if has_M
            Bc_k(:,:) = Ts * solveFactoredLinearSystemLU(L, U, Bc_k);
        else
            Bc_k(:,:) = Ts * Bc_k;
        end
    end
    if has_dist && (first_stage || k <= len_Bdc || factor_changed)
        Bdc_k(:,:) = Bdc(:,:,min(k,len_Bdc));
        if has_M
            Bdc_k(:,:) = Ts * solveFactoredLinearSystemLU(L, U, Bdc_k);
        else
            Bdc_k(:,:) = Ts * Bdc_k;
        end
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
