function mpc = update_mpc_general_lin_ineq_cnstr(mpc,Ci,Di,Ddi,...
                yi_min,yi_max)
arguments
    mpc
    Ci
    Di
    Ddi
    yi_min = []
    yi_max = []
end

update_grads = 0;

if ~isempty(Ci)
    mpc.Ci = Ci;
    update_grads = 1;
end

if ~isempty(Di)
    mpc.Di = Di;
    update_grads = 1;
end

if ~isempty(Ddi)   
    mpc.Ddi = Ddi;
end

if ~isempty(yi_min)   
    mpc.yi_min = yi_min;
end

if ~isempty(yi_max)    
    mpc.yi_max = yi_max;
end

% General Inequalites box constraints
if update_grads
    
    [mpc.gradYimin,mpc.gradYimax] = genGradY(mpc.Ci,mpc.Di,mpc.N,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.Nyi,mpc.nx,mpc.nu,mpc.nyi);
    
    if ~isempty(mpc.yi_min)
        [mpc.hessYimin,~] = genHessIneq(mpc.gradYimin);
    end
    if ~isempty(mpc.yi_max)
        [mpc.hessYimax,~] = genHessIneq(mpc.gradYimax);
    end

end

end