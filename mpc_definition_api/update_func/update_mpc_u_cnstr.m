%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_u_cnstr(mpc,u_min,u_max)
%
% Modifies the constraints limits on the control action
%
% In:
%   - mpc: CHRONOS mpc structure
%   - u_min (optional): nu column vector, lower bound constraint values on
%   the control action
%   - u_max (optional): nu column vector, upper bound constraint values on
%   the control action 
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_u_cnstr(mpc,u_min,u_max)

if ~isempty(u_min)    
    mpc.u_min = u_min;
end

if ~isempty(u_max)    
    mpc.u_max = u_max;
end

end