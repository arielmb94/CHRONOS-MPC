function [mpc] = update_mpc_ter_ingredients(mpc,P)

mpc.recompute_cost_hess = 1;

mpc.P = P;

mpc.hessTerminalCost(end-mpc.nx+1: end,end-mpc.nx+1 : end) = 2*P;

end