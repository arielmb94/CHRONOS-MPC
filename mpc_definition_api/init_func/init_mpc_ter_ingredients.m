function mpc = init_mpc_ter_ingredients(mpc,P,xN_ref_is_y)
% INIT_MPC_TER_INGREDIENTS Add a terminal-state cost using a supplied P.
%
%   mpc = INIT_MPC_TER_INGREDIENTS(mpc, P) adds the terminal penalty
%
%       J_terminal = (xN_ref - s_N)' * P * (xN_ref - s_N)
%
%   using the supplied terminal weight. P must be one nx-by-nx matrix; a
%   scalar is expanded as P*eye(nx). Zero P is stored as a zero weight;
%   P=[] is stored as empty. Both still enable terminal ingredients and
%   reset the terminal reference to zero; neither is a no-op. This function
%   leaves any existing rollout gain mpc.K unchanged.
%
%   mpc = INIT_MPC_TER_INGREDIENTS(mpc, P, true) enables the tracking-reference
%   shortcut. When xN_ref_in is empty, CHRONOS reuses the tracking reference
%   as the terminal-state reference. Use this option only when the tracking
%   output is the full state, y = s, with the same component order. An
%   explicitly supplied xN_ref_in always takes precedence.
%
%   By default, xN_ref_is_y is false and the user must pass the full-state
%   terminal reference separately as xN_ref_in to MPC_SOLVE.
%
%   Call this function after INIT_MPC_DYNAMICS and before BUILD_CHRONOS_MPC.
%
%   Inputs:
%     mpc          - CHRONOS MPC structure.
%     P            - Terminal weight, size nx-by-nx. A scalar scales eye(nx).
%                    Symmetric positive definiteness is assumed, not
%                    validated; empty and zero weights are not rejected.
%     xN_ref_is_y  - Optional flag. Set true to reuse the tracking
%                    reference when xN_ref_in is empty. Default: false.
%
%   Output:
%     mpc          - Updated CHRONOS MPC structure containing mpc.P, enabled
%                    terminal ingredients, and a zero terminal reference.
%                    The rollout gain mpc.K is unchanged.
%
%   Example - supply a terminal weight and an explicit runtime reference:
%
%       P = diag([10, 1]);
%       mpc = init_mpc_ter_ingredients(mpc, P);
%       % Pass xN_ref as xN_ref_in when calling MPC_SOLVE.

arguments
    mpc
    P
    xN_ref_is_y = 0
end

validate_matrix(P, mpc.nx, mpc.nx, 'P', true);

if size(P,3) > 1
    error('CHRONOS:DimensionMismatch', ...
        'Input "P" for terminal ingredients must contain one matrix page; supplied %d pages.', ...
        size(P,3));
end

if isscalar(P), P = P * eye(mpc.nx); end

mpc.ter_ingredients = 1;
if xN_ref_is_y
    mpc.xN_ref_is_y = 1;
else
    mpc.xN_ref_is_y = [];
end

mpc.xN_ref = zeros(mpc.nx,1);
mpc.P = P;

end
