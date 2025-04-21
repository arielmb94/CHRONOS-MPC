function mpc = update_mpc_DiffControl_cost(mpc,Rdu)

mpc.Rdu = Rdu;

[mpc.gradDiffCtlrRdu,mpc.hessDiffCtrlTerm] = genDiffControlGradHess(Rdu,...
    mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

mpc.recompute_cost_hess = 1;

end