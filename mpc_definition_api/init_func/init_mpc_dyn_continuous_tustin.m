function mpc = init_mpc_dyn_continuous_tustin(mpc, Ac, Bc, Bdc, Ts, M)
% INIT_MPC_DYN_CONTINUOUS_TUSTIN Define dynamics with the Tustin method.
%
%   mpc = INIT_MPC_DYN_CONTINUOUS_TUSTIN(mpc, Ac, Bc, Bdc, Ts) discretizes
%
%       x_dot = Ac*x + Bc*u + Bdc*d
%
%   Use [] for Bdc when there is no fixed known input d. An entirely zero
%   Bdc also omits the known-input term.
%
%   mpc = INIT_MPC_DYN_CONTINUOUS_TUSTIN(mpc, Ac, Bc, Bdc, Ts, M) uses
%   an optional mass matrix in the continuous model
%
%       M*x_dot = Ac*x + Bc*u + Bdc*d.
%
%   Pass Ac, Bc, and Bdc directly; this initializer incorporates M during
%   discretization, so no transformation by M is needed beforehand.
%   Omitted, empty, or scalar-zero M means I.
%
%   The Tustin model, with M = I by default, is
%
%       (M - (Ts/2)*Ac)*x_(k+1) = (M + (Ts/2)*Ac)*x_k
%                                + Ts*Bc*u_k + Ts*Bdc*d_k.
%
%   Ac, Bc, Bdc, and M may independently be constant matrices or arrays
%   with L pages. Each short input repeats its last page through horizon N;
%   pages beyond N are ignored.
%
%   Ts must be a positive scalar. Nondefault M must be invertible.
%
%   Call this function after INIT_CHRONOS_MPC and before adding costs or
%   constraints. INIT_MPC_DYNAMICS installs the full-state tracking output
%   y_k = x_k. Call INIT_MPC_OUTPUT afterward for a different output.
%
%   Inputs:
%     mpc     - CHRONOS MPC structure created by INIT_CHRONOS_MPC.
%     Ac      - State coefficient, size nx-by-nx or nx-by-nx-by-L.
%     Bc      - Control coefficient, size nx-by-nu or nx-by-nu-by-L.
%     Bdc     - Fixed-known-input coefficient, size nx-by-nd or
%               nx-by-nd-by-L. Pass [] to omit the term.
%     Ts      - Positive scalar sampling interval.
%     M       - Optional mass matrix, size nx-by-nx or nx-by-nx-by-L.
%               Default: I for omitted, empty, or scalar-zero M.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - discretize a model without a known input:
%
%       mpc = init_chronos_mpc(20);
%       mpc = init_mpc_dyn_continuous_tustin(mpc, Ac, Bc, [], Ts);
%
%   Example - combine staged state and mass matrices with constant Bc:
%
%       Ac = cat(3, Ac_1, Ac_2);
%       M = cat(3, M_1, M_2, M_3);
%       mpc = init_mpc_dyn_continuous_tustin(mpc, Ac, Bc, [], Ts, M);

arguments
    mpc
    Ac
    Bc
    Bdc
    Ts
    M = []
end

% Normalize omitted optional terms before validating their dimensions.
if ~isempty(Bdc) && ~any(Bdc(:))
    Bdc = [];
end
if isempty(M) || (isscalar(M) && M == 0)
    M = [];
end
has_dist = ~isempty(Bdc);
has_M = ~isempty(M);

% Validate the original coefficients before discretization arithmetic.
nx = size(Ac,1);
validate_matrix(Ac, nx, nx, 'Ac');
nu = size(Bc,2);
validate_matrix(Bc, nx, nu, 'Bc');

len_Bdc = 0;
if has_dist
    nd = size(Bdc,2);
    validate_matrix(Bdc, nx, nd, 'Bdc');
    len_Bdc = size(Bdc,3);
end
len_M = 0;
if has_M
    validate_matrix(M, nx, nx, 'M');
    len_M = size(M,3);
end

% Discretize supplied stages; INIT_MPC_DYNAMICS fills the remaining horizon.
len_Ac = size(Ac,3);
len_Bc = size(Bc,3);
num_stages = min(mpc.N, max([len_Ac, len_Bc, len_Bdc, len_M]));
Ad = zeros(nx, nx, num_stages);
Bd = zeros(nx, nu, num_stages);
Bdd = [];
if has_dist
    Bdd = zeros(nx, nd, num_stages);
end
I = eye(nx);
Ac_k = Ad(:,:,1);
Bc_k = Bd(:,:,1);
if has_dist
    Bdc_k = Bdd(:,:,1);
end
M_k = I;
W = Ad(:,:,1);
rhs_A = Ad(:,:,1);
L = Ad(:,:,1);
U = Ad(:,:,1);

for k = 1:num_stages
    first_stage = (k == 1);
    factor_changed = first_stage || k <= len_Ac || (has_M && k <= len_M);
    if factor_changed
        Ac_k(:,:) = Ac(:,:,min(k,len_Ac));
        if has_M
            M_k(:,:) = M(:,:,min(k,len_M));
        end
        W(:,:) = M_k - (Ts/2) * Ac_k;
        rhs_A(:,:) = M_k + (Ts/2) * Ac_k;
        [L(:,:), U(:,:)] = factorizeLinearSystemLU(W);
        Ac_k(:,:) = solveFactoredLinearSystemLU(L, U, rhs_A);
    end
    if factor_changed || k <= len_Bc
        Bc_k(:,:) = Bc(:,:,min(k,len_Bc));
        Bc_k(:,:) = solveFactoredLinearSystemLU(L, U, Ts * Bc_k);
    end
    if has_dist && (factor_changed || k <= len_Bdc)
        Bdc_k(:,:) = Bdc(:,:,min(k,len_Bdc));
        Bdc_k(:,:) = solveFactoredLinearSystemLU(L, U, Ts * Bdc_k);
    end

    Ad(:,:,k) = Ac_k;
    Bd(:,:,k) = Bc_k;
    if has_dist
        Bdd(:,:,k) = Bdc_k;
    end
end

mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end
