% INIT_MPC_TER_INGREDIENTS_DLQR Add a DLQR terminal-state cost.
%
%   mpc = INIT_MPC_TER_INGREDIENTS_DLQR(mpc, Qx, Ru) computes
%
%       [K, P] = dlqr(A_N, B_N, Qx, Ru)
%
%   from the final horizon dynamics and adds the terminal penalty
%
%       J_terminal = (xN_ref - s_N)' * P * (xN_ref - s_N)
%
%   Store representative dynamics with INIT_MPC_DYNAMICS before calling
%   this function. The calculation uses terminal mpc.A(:,:,mpc.N) and
%   mpc.B(:,:,mpc.N).
%
%   mpc = INIT_MPC_TER_INGREDIENTS_DLQR(mpc, Qx, Ru, true) enables the
%   tracking-reference shortcut. When xN_ref_in is empty, CHRONOS reuses
%   the tracking reference as the terminal-state reference. Use this option
%   only when the tracking output is the full state, y = s, with the same
%   component order. An explicitly supplied xN_ref_in always takes
%   precedence.
%
%   By default, xN_ref_is_y is false and the user must pass the full-state
%   terminal reference separately as xN_ref_in to MPC_SOLVE.
%
%   Call this function after INIT_MPC_DYNAMICS and before
%   BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc          - CHRONOS MPC structure.
%     Qx           - State weight passed to DLQR, size nx-by-nx.
%     Ru           - Control weight passed to DLQR, size nu-by-nu.
%     xN_ref_is_y  - Optional flag. Set true to reuse the tracking
%                    reference when xN_ref_in is empty. Default: false.
%
%   Output:
%     mpc          - Updated CHRONOS MPC structure containing mpc.K and
%                    mpc.P.
%
%   Example - require an explicit terminal-state reference at runtime:
%
%       Qx = eye(nx);
%       Ru = 0.1 * eye(nu);
%       mpc = init_mpc_ter_ingredients_dlqr(mpc, Qx, Ru);
%       % Pass xN_ref as xN_ref_in when calling MPC_SOLVE.
function mpc = init_mpc_ter_ingredients_dlqr(mpc,Qx,Ru,...
                                             xN_ref_is_y)
arguments
    mpc
    Qx
    Ru
    xN_ref_is_y = 0
end

if isempty(Qx) || ~any(Qx(:)) || isempty(Ru) || ~any(Ru(:))
    return;
end

validate_matrix(Qx, mpc.nx, mpc.nx, 'Qx', true);
validate_matrix(Ru, mpc.nu, mpc.nu, 'Ru', true);

if size(Qx,3) > 1
    error('CHRONOS:DimensionMismatch', ...
        'Input "Qx" for DLQR must contain one matrix page; supplied %d pages.', ...
        size(Qx,3));
end
if size(Ru,3) > 1
    error('CHRONOS:DimensionMismatch', ...
        'Input "Ru" for DLQR must contain one matrix page; supplied %d pages.', ...
        size(Ru,3));
end

if isscalar(Qx), Qx = Qx * eye(mpc.nx); end
if isscalar(Ru), Ru = Ru * eye(mpc.nu); end

mpc.ter_ingredients = 1;
if xN_ref_is_y
    mpc.xN_ref_is_y = 1;
else
    mpc.xN_ref_is_y = [];
end

[K,P] = dlqr(mpc.A(:,:,mpc.N),mpc.B(:,:,mpc.N),Qx,Ru);

mpc.xN_ref = zeros(mpc.nx,1);

mpc.K = K;
mpc.P = P;

end
