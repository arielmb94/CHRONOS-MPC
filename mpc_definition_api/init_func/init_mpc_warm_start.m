%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,d,x_ref,dh,x0)
%
% Computes a feasible value for the optimization variables vector x0.
%
% Example Use:
%
%   - mpc problem has no disturbance inputs, terminal set constraints
%   are not enabled:
%           x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev) 
%   - mpc problem has no disturbance inputs, terminal set constraints
%   are enabled:
%           x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,[],x_ref) 
%   - mpc problem has no disturbance inputs, terminal set constraints
%   are not enabled, initial feasible guess is provided:
%           x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,[],[],[],x0)
%   - mpc problem has disturbance inputs on user defined constrained signal 
%   h but not on system dynamics, terminal set constraints are not enabled
%           x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,[],[],dh)
%
% In:
%   - mpc: CHRONOS mpc structure
%   - s_prev: nx column vector, initial value for the state vector
%   - u_prev: nu column vector, initial value for the control action vector
%   - d (optional): nd column vector or Nd column vector, vector of known 
%   disturbance inputs, required when the system dynamics have a defined Bd
%   or Dd matrices
%   - x_ref (optional): nx column vector, reference value for the terminal 
%   ingredients set constraint, only required when the terminal set 
%   constraint is enabled
%   - dh (optional): ndh column vector, vector of known disturbance inputs 
%   for the custom user defined signal h, required if the user has defined
%   custom constraints with disturbance input on the signal h. 
%   - x0 (optional): Nx + Nu column vector, initial guess for a feasible 
%   value for the optimization variables vector x0. If the x0 argument is
%   defined, init_mpc_warm_start() will try to find the closest feasible
%   point to x0
%
% Out:
%
%   - x_mpc: computed feasible value for the optimization variables vector
%   x0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,d,x_ref,dh,x0)
arguments
mpc
s_prev
u_prev
d = [];
x_ref = [];
dh = [];
x0 = [];
end

% enable warm_starting boolean
mpc.warm_starting = 1;

if isempty(x0)
    % initialize optimization vector to 0
    x0 = zeros(mpc.Nx+mpc.Nu,1);
end

% update b matrix from equality condition
mpc = update_mpc_beq(mpc,s_prev,d);

% find feasible point
[x_mpc,iter] = feas_solve(x0,mpc,s_prev,u_prev,d,x_ref,dh);

% disable warm_starting boolean
mpc.warm_starting = 0;

end