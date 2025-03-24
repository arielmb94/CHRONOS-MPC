function x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,d,x_ref,di,x0)
arguments
mpc
s_prev
u_prev
d = [];
x_ref = [];
di = [];
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
[x_mpc,iter] = feas_solve(x0,mpc,s_prev,u_prev,d,x_ref,di);

% disable warm_starting boolean
mpc.warm_starting = 0;

end