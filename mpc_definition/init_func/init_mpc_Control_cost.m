function mpc = init_mpc_Control_cost(mpc,Ru,ru)
arguments
    mpc
    Ru = []
    ru = [];
end

mpc.Ru = Ru;
mpc.ru = ru;

% Quadratic Cost
if ~isempty(Ru)

    if isempty(mpc.hessCost)
        mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
    end

    [mpc.gradCtlrRu,mpc.hessCtrlTerm] = genControlGradHess(Ru,mpc.N_ctr_hor,...
        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

    mpc.hessCost = mpc.hessCost + mpc.hessCtrlTerm;
end

% Linear Cost
if ~isempty(ru)

    mpc.gradCtlrru = genControlLPGrad(ru,mpc.N_ctr_hor,mpc.Nx,mpc.Nu,...
                    mpc.nx,mpc.nu);
    
end

end