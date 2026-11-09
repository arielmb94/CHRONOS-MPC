function mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, Bdc, Ts, M)
% INIT_MPC_DYN_CONTINUOUS_FORWARD Define dynamics with forward Euler.
%
%   mpc = INIT_MPC_DYN_CONTINUOUS_FORWARD(mpc, Ac, Bc, Bdc, Ts) discretizes
%
%       x_dot = Ac*x + Bc*u + Bdc*d
%
%   Use [] for Bdc when there is no fixed known input d. An entirely zero
%   Bdc also omits the known-input term.
%
%   mpc = INIT_MPC_DYN_CONTINUOUS_FORWARD(mpc, Ac, Bc, Bdc, Ts, M) uses
%   an optional mass matrix in the continuous model
%
%       M*x_dot = Ac*x + Bc*u + Bdc*d.
%
%   Pass Ac, Bc, and Bdc directly; this initializer incorporates M during
%   discretization, so no transformation by M is needed beforehand.
%   Omitted, empty, or scalar-zero M means I.
%
%   The forward Euler model, with M = I by default, is
%
%       x_(k+1) = Ad_k*x_k + Bd_k*u_k + Bdd_k*d_k
%       Ad = I + Ts*M^(-1)*Ac, Bd = Ts*M^(-1)*Bc,
%       Bdd = Ts*M^(-1)*Bdc.
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
%       mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, [], Ts);
%
%   Example - combine staged state and mass matrices with constant Bc:
%
%       Ac = cat(3, Ac_1, Ac_2);
%       M = cat(3, M_1, M_2, M_3);
%       mpc = init_mpc_dyn_continuous_forward(mpc, Ac, Bc, [], Ts, M);

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
L = Ad(:,:,1);
U = Ad(:,:,1);

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

    Ad(:,:,k) = Ac_k;
    Bd(:,:,k) = Bc_k;
    if has_dist
        Bdd(:,:,k) = Bdc_k;
    end
end

mpc = init_mpc_dynamics(mpc, Ad, Bd, Bdd);
end
