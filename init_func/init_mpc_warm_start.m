function x_mpc = init_mpc_warm_start(mpc,s_prev,u_prev,d,x_ref,di)
arguments
mpc
s_prev
u_prev
d = [];
x_ref = [];
di = [];
end

% enable warm_starting boolean
mpc.warm_starting = 1;

% initialize optimization vector to 0
x0 = zeros(mpc.Nx+mpc.Nu,1);

% find feasible point
x_mpc = feas_solve(x0,mpc,s_prev,u_prev,d,x_ref,di);

% disable warm_starting boolean
mpc.warm_starting = 0;

end