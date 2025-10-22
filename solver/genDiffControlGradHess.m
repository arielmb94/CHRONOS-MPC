function mpc = genDiffControlGradHess(mpc,R,N,nx,nu)

for k = 0:N-1

    switch k
        case 0
            % Write in uo
            mpc.gradDiffCtlrR(1:nu,1:nu) = R;
            mpc.gradDiffCtlr(1:nu,1:nu) = eye(nu);
        otherwise
            % write in uk-1
            mpc.gradDiffCtlrR(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -R;
            mpc.gradDiffCtlr(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -eye(nu);
            % write in uk
            mpc.gradDiffCtlrR(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                R;
            mpc.gradDiffCtlr(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                eye(nu);

    end
end

%Grad term is then DeltaJu = gradCtlrR*deltaU

mpc.hessDiffCtrlTerm(:,:) = mpc.gradDiffCtlrR*mpc.gradDiffCtlr';


end