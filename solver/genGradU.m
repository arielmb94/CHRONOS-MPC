function [gradUmin,gradUmax] = genGradU(N,Nx,Nu,nx,nu)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% u_min - uk <= 0 
% uk - umax <= 0 

gradUmax = zeros(Nx+Nu,Nu);

for k = 0:N-1
    switch k
        case 0
            gradUmax(1:nu,1:nu) = eye(nu);
        otherwise
            gradUmax(nu*k+nx*k+1: nu*k+nx*k+nu, nu*(k)+1:nu*(k+1)) = eye(nu);
    end
    
end

gradUmin = -gradUmax;
end