function [gradDeltaUmin,gradDeltaUmax] = genGradDeltaU(N,Nx,Nu,nx,nu)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% u_min - uk + uk-1 <= 0 
% uk - uk-1 - umax <= 0 

gradDeltaUmax = zeros(Nx+Nu,Nu);

for k = 0:N-1

    switch k
        case 0
            % Write in uo
            gradDeltaUmax(1:nu,1:nu) = eye(nu);
        otherwise
            % write in uk-1
            gradDeltaUmax(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1),...
                nu*(k)+1:nu*(k+1)) = -eye(nu);
            % write in uk
            gradDeltaUmax(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k),...
                nu*(k)+1:nu*(k+1)) = eye(nu);

    end
end

gradDeltaUmin = -gradDeltaUmax;
end