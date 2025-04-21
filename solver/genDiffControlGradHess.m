function [gradDiffCtlrR,hessDiffCtrlTerm] = genDiffControlGradHess(R,N,Nx,Nu,nx,nu)


gradDiffCtlrR = zeros(Nx+Nu,Nu);
gradDiffCtlr = zeros(Nx+Nu,Nu);
for k = 0:N-1

    switch k
        case 0
            % Write in uo
            gradDiffCtlrR(1:nu,1:nu) = R;
            gradDiffCtlr(1:nu,1:nu) = eye(nu);
        otherwise
            % write in uk-1
            gradDiffCtlrR(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -R;
            gradDiffCtlr(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -eye(nu);
            % write in uk
            gradDiffCtlrR(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                R;
            gradDiffCtlr(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                eye(nu);

    end
end

%Grad term is then DeltaJu = gradCtlrR*deltaU

hessDiffCtrlTerm = gradDiffCtlrR*gradDiffCtlr';


end