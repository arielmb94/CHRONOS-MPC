%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_state_cnstr(mpc,x_min,x_max)
%
% Modifies the constraints limits on the state vector x.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - x_min (optional): nx column vector, lower bound constraint values on 
%   the state vector
%   - x_max (optional): nx column vector, upper bound constraint values on 
%   the state vector
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_state_cnstr(mpc,x_min,x_max)

if ~isempty(x_min)    
    mpc.x_min = x_min;
end

if ~isempty(x_max)    
    mpc.x_max = x_max;
end

end