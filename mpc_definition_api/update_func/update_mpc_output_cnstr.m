%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   mpc = update_mpc_output_cnstr(mpc,y_min,y_max)
%
% Modifies the constraints limits on the tracking signal y.
%
% In:
%   - mpc: CHRONOS mpc structure
%   - y_min (optional): ny column vector, lower bound constraint values on 
%   the tracking signal
%   - y_max (optional): ny column vector, upper bound constraint values on 
%   the tracking signal
%
% Out:
%   - mpc: updated CHRONOS mpc structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = update_mpc_output_cnstr(mpc,y_min,y_max)

if ~isempty(y_min)    
    mpc.y_min = y_min;
end

if ~isempty(y_max)    
    mpc.y_max = y_max;
end

end