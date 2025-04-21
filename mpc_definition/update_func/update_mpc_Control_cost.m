function mpc = update_mpc_Control_cost(mpc,Ru,ru)
arguments
    mpc
    Ru = []
    ru = []
end

if ~isempty(Ru)
    
    mpc.Ru = Ru;
    
    [mpc.gradCtlrRu,mpc.hessCtrlTerm] = genControlGradHess(Ru,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
    
    mpc.recompute_cost_hess = 1;

end

if ~isempty(ru)

    mpc.ru = ru;

    mpc.gradCtlrru = genControlLPGrad(ru,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,...
                    mpc.nx,mpc.nu);

end

end