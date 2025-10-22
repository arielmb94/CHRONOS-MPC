%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max)
%
% Modifies the box constraints limits on the terminal state vector xN.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - x_ter_min (optional): nx column vector, lower bound constraint values
%   on the terminal state vector
%   - x_ter_max (optional): nx column vector, upper bound constraint values
%   on the terminal state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_ter_state_cnstr(mpc,x_ter_min,x_ter_max)

if ~isempty(x_ter_min)    
    mpc.s_ter_cnstr.min = x_ter_min;
end

if ~isempty(x_ter_max)    
    mpc.s_ter_cnstr.max = x_ter_max;
end

end