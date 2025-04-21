function grad_J = grad_f0_MPC(mpc,err,deltaU,U,grad_ter,z)

    grad_J = zeros(mpc.Nx+mpc.Nu,1);

    if ~isempty(mpc.gradErrQe)
        grad_J = grad_J - mpc.gradErrQe*err;
    end    

    if ~isempty(mpc.gradDiffCtlrRdu)
        grad_J = grad_J + mpc.gradDiffCtlrRdu*deltaU;
    end  

    if ~isempty(mpc.gradCtlrRu)
        grad_J = grad_J + mpc.gradCtlrRu*U;
    end

    if ~isempty(mpc.gradCtlrru)
        grad_J = grad_J + mpc.gradCtlrru;
    end

    if mpc.ter_ingredients
        grad_J(end-mpc.nx+1:end) = grad_J(end-mpc.nx+1:end) - grad_ter;
    end  

    if ~isempty(mpc.gradPerfQz)
        grad_J = grad_J + mpc.gradPerfQz*z;
    end  

    if ~isempty(mpc.gradPerfqz)
        grad_J = grad_J + mpc.gradPerfqz;
    end

end