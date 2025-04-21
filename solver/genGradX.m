function [gradXmin,gradXmax] = genGradX(N,N_h_ctr,Nx,Nu,nx,nu)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% x_min - xk <= 0 
% xk - xmax <= 0 

gradXmax = zeros(Nx+Nu,Nx);
for k = 1:N
    if k < N_h_ctr
        gradXmax(nu + nx*(k-1)+ nu*(k-1) + 1 : nu + nx*(k)+ nu*(k-1),...
            nx*(k-1)+1:nx*(k)) = eye(nx);

    else
        gradXmax(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
            nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
            nx*(k-1)+1:nx*(k)) = eye(nx);

    end
end

gradXmin = -gradXmax;
end