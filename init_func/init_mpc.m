function [mpc] = init_mpc(N,N_ctr_hor)
arguments
    N = 2
    N_ctr_hor = 0
end

mpc.N = N;
if N_ctr_hor && N_ctr_hor > N
    mpc.N_ctr_hor = N;
elseif N_ctr_hor
    mpc.N_ctr_hor = N_ctr_hor;
else
    mpc.N_ctr_hor = N;
end
mpc.Qe = [];
mpc.Rdu = [];
mpc.Ru = [];
mpc.ru = [];
mpc.A = [];
mpc.B = [];
mpc.Bd = [];
mpc.C = [];
mpc.D = [];
mpc.Dd = [];
mpc.Qz = [];
mpc.qz = [];
mpc.Cz = [];
mpc.Dz = [];
mpc.Ddz = [];
mpc.Ci = [];
mpc.Di = [];
mpc.Ddi = [];
mpc.nx = 0;
mpc.nu = 0;
mpc.nd = 0;
mpc.ny = 0;
mpc.ndz = 0;
mpc.nz = 0;
mpc.ndi = 0;
mpc.nyi = 0;
mpc.Nx = 0;
mpc.Nu = 0;
mpc.Nd = 0;
mpc.Ny = 0;
mpc.Nz = 0;
mpc.Nyi = 0;
mpc.Ndz = 0;
mpc.Ndi = 0;
mpc.x_min = [];
mpc.x_max = [];
mpc.x_ter_min = [];
mpc.x_ter_max = [];
mpc.u_min = [];
mpc.u_max = [];
mpc.du_min = [];
mpc.du_max = [];
mpc.y_min = [];
mpc.y_max = [];
mpc.yi_min = [];
mpc.yi_max = [];
mpc.Aeq = [];
mpc.beq = [];
mpc.hessCost = [];
mpc.gradErrQe = [];
mpc.hessErrTerm = [];
mpc.gradDiffCtlrRdu = [];
mpc.hessDiffCtrlTerm = [];
mpc.gradCtlrRu = [];
mpc.gradCtlrru = [];
mpc.hessCtrlTerm = [];
mpc.gradPerfQz = [];
mpc.hessPerfTerm = [];
mpc.gradXmin = [];
mpc.gradXmax = [];
mpc.hessXmin = [];
mpc.hessXmax = [];
mpc.gradUmin = [];
mpc.gradUmax = [];
mpc.hessUmin = [];
mpc.hessUmax = [];
mpc.gradDeltaUmin = [];
mpc.gradDeltaUmax = [];
mpc.hessDeltaUmin = [];
mpc.hessDeltaUmax = [];
mpc.m = 0;
mpc.t = 50;
mpc.Beta = 0.75;
mpc.min_l = 0.99;
mpc.eps = 1e-4;
mpc.ter_ingredients = 0;
mpc.ter_constraint = 0;
mpc.x_ref_is_y = 0;
mpc.P = [];
mpc.hessTerminalCost = [];
mpc.recompute_cost_hess = 0;

end