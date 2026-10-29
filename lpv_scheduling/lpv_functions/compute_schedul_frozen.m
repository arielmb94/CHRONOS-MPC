%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Pk_frozen = compute_schedul_frozen(mpc, x_prev, sched_fun)
%
% Computes the scheduling trajectory by evaluating the current measured 
% state and keeping them constant (frozen) along the entire prediction 
% horizon.
%
% In:
%   - mpc: CHRONOS mpc structure.
%   - x_prev: nx column vector, current measured or estimated system state.
%   - sched_fun: function handle (@(x)), returns n_rho scheduling parameters 
%   evaluated at state x.
%
% Out:
%   - Pk_frozen: (N * n_rho) column vector, frozen scheduling trajectory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pk_frozen = compute_schedul_frozen(mpc, x_prev, sched_fun)
    % Compute scheduling variables for current instant only
    rho_current = sched_fun(x_prev);
    
    % Freeze (repeat) vector along prediction horizon
    Pk_frozen = repmat(rho_current, mpc.N, 1);
end