% Computes gradient for linear cost terms of the form:
% qz^t*z
% where y is a linear output equation:
% z = Cz*s + Dz*u + Ddz*d
%
% with:
%
% grad = qz

function gradq = genGenPerfLPGrad(q,C,D,N,N_h_ctr,Nx,Nu,nx,nu)

% gradQ = Grad(y)*Q
gradq = zeros(Nx+Nu,1);

% Precompute vector multiplications
Dq = D'*q;
Cq = C'*q;
CDq = [C';D']*q;

for k = 0:N-1

    if k == 0

        if D==0
            % Do nothing
        else

            gradq(1:nu) = Dq;
            
        end

    elseif k < N_h_ctr

        if D==0 % -> Ny is (N-1)*ny

            gradq(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k-1))...
                    = Cq;
  
        else    % -> Ny is N*ny

            gradq(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k)) ... 
                    = CDq;

        end

    else

        if D==0

            gradq(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx) = Cq;
       
        else

            gradq(nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + 1  : ...
                nu + (nx+nu)*(N_h_ctr-1) + nx*(k-N_h_ctr) + nx) = Cq;
            
            % gradYmax = D' -> u^N_h_ctr-1
            gradq(nu + (nx+nu)*(N_h_ctr-2) + nx + 1 : ...
                nu + (nx+nu)*(N_h_ctr-1)) =  Dq;

        end
    end   
end

end