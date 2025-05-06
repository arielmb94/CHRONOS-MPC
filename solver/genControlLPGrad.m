function gradCtlrru = genControlLPGrad(r,N,Nx,Nu,nx,nu)

gradCtlrru = zeros(Nx+Nu,1);
for k = 0:N-1

    switch k
        case 0
            % Write in uo
            gradCtlrru(1:nu) = r;
        otherwise
            % write in uk
            gradCtlrru(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k)) = r;
    end
    
end

end