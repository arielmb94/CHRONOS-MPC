% Computes output vectors of the form: y = Cs + Du + Dd
function y = get_lin_out(s,u,d,nx,nu,ny,nd,N,N_h_ctr,Ny,C,D,Dd,Nd)

    if isempty(d) || isempty(Dd) % No disturbance term

        y = zeros(Ny,1);
        for k = 1:Ny/ny
            
            if D == 0 % -> Ny is (N-1)*ny

                % ignore k = 0, start from k = 1
                sk = s(((k-1)*nx+1:k*nx)+nx);
                
                y((k-1)*ny+1:k*ny) = C*sk;

            else      % -> Ny is (N)*ny

                sk = s((k-1)*nx+1:k*nx);
                if k < N_h_ctr
                    uk = u((k-1)*nu+1:k*nu);
                else
                    uk = u(end-nu+1:end);
                end

                y((k-1)*ny+1:k*ny) = C*sk+D*uk;

            end
        end
    
    else % There are disturbance terms

        % copy d N times if only given for a sampling instance
        if length(d)<Nd
            d = repmat(d,N,1);
        end

        y = zeros(Ny,1);
        for k = 1:Ny/ny
            
            if D == 0 % -> Ny is (N-1)*ny

                % ignore k = 0, start from k = 1
                sk = s(((k-1)*nx+1:k*nx)+nx);
                dk = d(((k-1)*nd+1:k*nd)+nd);
                
                y((k-1)*ny+1:k*ny) = C*sk+Dd*dk;

            else      % -> Ny is (N)*ny

                sk = s((k-1)*nx+1:k*nx);
                if k < N_h_ctr
                    uk = u((k-1)*nu+1:k*nu);
                else
                    uk = u(end-nu+1:end);
                end
                dk = d((k-1)*nd+1:k*nd);

                y((k-1)*ny+1:k*ny) = C*sk+D*uk+Dd*dk;

            end
        end
        
    end
end