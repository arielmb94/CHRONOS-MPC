function mpc = update_mpc_LinPerf_cost(mpc,Cz,Dz,Ddz,Qz,qz)

updatelp = 0;
updateqp = 0;

if ~isempty(Cz)
    mpc.Cz = Cz;

    updatelp = 1;
    updateqp = 1;
end

if ~isempty(Dz)
    mpc.Dz = Dz;

    updatelp = 1;
    updateqp = 1;
end

if ~isempty(Ddz)   
    mpc.Ddz = Ddz;
end

if ~isempty(Qz)   
    mpc.Qz = Qz;

    updateqp = 1;
end

if ~isempty(qz)   
    mpc.qz = qz;

    updatelp = 1;
end

% Update Quatric cost term gradient and Hessian
if ~isempty(mpc.Qz) && updateqp

    mpc.recompute_cost_hess = 1;

    [mpc.gradPerfQz,mpc.hessPerfTerm] = genLinOutGradHess(mpc.Qz,mpc.Cz,mpc.Dz,...
        mpc.N,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,mpc.Nz,mpc.nx,mpc.nu,mpc.nz);

end

% Update Linear cost term gradient
if ~isempty(mpc.qz) && updatelp

    mpc.gradPerfqz = genGenPerfLPGrad(mpc.qz,mpc.Cz,mpc.Dz,mpc.N,mpc.N_ctr_hor,...
                        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

end
    
end