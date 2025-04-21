function [gradYmin,gradYmax] = genGradY(C,D,N,N_h_ctr,Nx,Nu,Ny,nx,nu,ny)
%optimization variables are stored as [u0 x1 u1 x2 ... xn-1 un-1 xn]^T

gradYmax = zeros(Nx+Nu,Ny);
for k = 0:N-1

    if k == 0

        if D == 0
            % Do Nothing
        else
            % gradYmax = D'
            gradYmax(1:nu,1:ny) = D';
        end

    elseif k < N_h_ctr

        if D == 0 % -> Ny is (N-1)*ny
            % gradYmax = C'
            gradYmax(nu + 1 + nx*(k-1) + nu*(k-1): nu + nx*(k)+ nu*(k-1),...
                (ny*(k)+1:ny*(k+1))-ny) =  C';
        else      % -> Ny is N*ny
            % gradYmax = [C';D']
            gradYmax(nu + 1 + nx*(k-1) + nu*(k-1): nu + nx*(k)+ nu*(k),...
                ny*(k)+1:ny*(k+1)) = [ C' ; D'];
        end

    else

        if D == 0 % -> Ny is (N-1)*ny
            % gradYmax = C'
            gradYmax(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))-ny) =  C';
        else      % -> Ny is N*ny
            % gradYmax = C' -> x^k
            gradYmax(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))) =  C';
            % gradYmax = D' -> u^N_h_ctr-1
            gradYmax(nu + (nx+nu)*(N_h_ctr-2) + nx + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) ,...
                (ny*(k)+1:ny*(k+1))) =  D';
        end
    end

end

gradYmin = -gradYmax;
end