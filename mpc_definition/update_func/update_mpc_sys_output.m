function mpc = update_mpc_sys_output(mpc,C,D,Dd,Qe,y_min,y_max)
arguments
    mpc
    C
    D
    Dd
    Qe = []
    y_min = []
    y_max = []
end

if ~isempty(C)
    mpc.C = C;
end

if ~isempty(D)
    mpc.D = D;
end

if ~isempty(Dd)   
    mpc.Dd = Dd;
end

if ~isempty(y_min)   
    mpc.y_min = y_min;
end

if ~isempty(y_max)    
    mpc.y_max = y_max;
end

% If tracking Cost exists, update gradients
if ~isempty(mpc.Qe)
    mpc = update_mpc_Tracking_cost(mpc,Qe);
end

% If Outputs constraint exists, update box constraints gradients
if ~isempty(mpc.y_min) ||  ~isempty(mpc.y_max)
    [mpc.gradYmin,mpc.gradYmax] = genGradY(mpc.C,mpc.D,mpc.N,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

    if ~isempty(mpc.y_min)
        [mpc.hessYmin,~] = genHessIneq(mpc.gradYmin);
    end
    if ~isempty(mpc.y_max)
        [mpc.hessYmax,~] = genHessIneq(mpc.gradYmax);
    end
end

end