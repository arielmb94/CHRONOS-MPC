% Computes gradient and hessian for cost terms of the form:
% 0.5*y'Qy
% where y is a linear output equation:
% y = Cs + Du + Dd
%
% with:
%
% grad = Grad(y)*Q*y
% hess = Grad(y)*Q*Delta_y'
%
% For error cost:
% e = r-y = r-(Cs + Du + Dd)
% grad = -[C'; D']*Q*e
% hess = [C' ; D']*Q*[C D]
%
% For general linear perf. cost:
% z = Cs + Du + Dd
% grad = [C'; D']*Q*z
% hess = [C' ; D']*Q*[C D]

function [gradQ,hess] = genLinOutGradHess(Q,C,D,N,N_h_ctr,Nx,Nu,Ny,nx,nu,ny)

% grad = Grad(y)*Q*y
% gradQ = Grad(y)*Q
grad = zeros(Nx+Nu,Ny);
gradQ = zeros(Nx+Nu,Ny);
for k = 0:N-1

    if k == 0

        if D==0
            % Do nothing
        else
            % gradQ = D'*Q
            grad(1:nu, 1:ny) = D';

            gradQ(1:nu, 1:ny) = grad(1:nu, 1:ny)*Q;
            
        end

    elseif k < N_h_ctr

        if D==0 % -> Ny is (N-1)*ny

            % gradQ = C'*Q
            grad(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k-1),...
                (ny*(k)+1:ny*(k+1))-ny) = C';

            gradQ(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k-1),...
                (ny*(k)+1:ny*(k+1))-ny) = grad(nu + 1 + nx*(k-1) + nu*(k-1) ...
                : nu + nx*(k)+ nu*(k-1),(ny*(k)+1:ny*(k+1))-ny)*Q;
  
        else    % -> Ny is N*ny

            % gradQ = [C';D']*Q
            grad(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                ny*(k)+1:ny*(k+1)) = [ C' ; D'];

            gradQ(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),...
                ny*(k)+1:ny*(k+1)) = grad(nu + 1 + nx*(k-1) + nu*(k-1)...
                : nu + nx*(k)+ nu*(k), ny*(k)+1:ny*(k+1))*Q;

        end

    else

        if D==0

            % gradQ = C'*Q
            grad(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))-ny) =  C';

            gradQ(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))-ny) = ...
                grad(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))-ny)*Q;
       
        else

            % gradYmax = C' -> x^k
            grad(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))) =  C';

            gradQ(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1))) = ...
                grad(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx,...
                (ny*(k)+1:ny*(k+1)))*Q;
            
            % gradYmax = D' -> u^N_h_ctr-1
            grad(nu + (nx+nu)*(N_h_ctr-2) + nx + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) ,...
                (ny*(k)+1:ny*(k+1))) =  D';

            gradQ(nu + (nx+nu)*(N_h_ctr-2) + nx + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) ,...
                (ny*(k)+1:ny*(k+1))) =  grad(nu + (nx+nu)*(N_h_ctr-2) + nx + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) , (ny*(k)+1:ny*(k+1)))*Q;

        end
    end   
end

hess = gradQ*grad';
end