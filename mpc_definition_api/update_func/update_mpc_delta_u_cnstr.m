%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_delta_u_cnstr(mpc,du_min,du_max)
%
% Modifies the constraint limits on the variation of the control action 
% between consecutive sampling instances.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - du_min (optional): nu column vector, lower bound constraint values on
%   the control variation between sampling instances
%   - du_max (optional): nu column vector, upper bound constraint values on
%   the control variation between sampling instances 
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_delta_u_cnstr(mpc,du_min,du_max)

if ~isempty(du_min)    
    mpc.du_min = du_min;
end

if ~isempty(du_max)    
    mpc.du_max = du_max;
end

end