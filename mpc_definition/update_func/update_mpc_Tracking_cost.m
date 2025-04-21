function mpc = update_mpc_Tracking_cost(mpc,Qe)

mpc.recompute_cost_hess = 1;

if ~isempty(Qe)   
    mpc.Qe = Qe;   
end

[mpc.gradErrQe,mpc.hessErrTerm] = genLinOutGradHess(mpc.Qe,mpc.C,mpc.D,...
            mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

end