function mpc = update_mpc_beq(mpc,x_prev,d,A0,Bd0,Pk,n_rho,varargin)
    if nargin < 4 || isempty(A0)
        A0 = mpc.A;
    end
    if nargin < 5 || isempty(Bd0)
        Bd0 = mpc.Bd;
    end
    
    if isempty(d) || isempty(Bd0)
        Ak = A0;
        
        if nargin > 5 && ~isempty(Pk)
            for i = 1:n_rho
                Ak = Ak + varargin{i} * Pk(i);
            end
        end
        
        mpc.beq(1:mpc.nx) = -Ak*x_prev;
    else
        for k = 0:mpc.N-1
            Ak = A0;
            Bdk = Bd0;
            
            if nargin > 5 && ~isempty(Pk)
                idx_base = k * n_rho;
                for i = 1:n_rho
                    Ak = Ak + varargin{i} * Pk(idx_base + i);
                    if length(varargin)/n_rho == 3
                        Bdk = Bdk + varargin{i + 2*n_rho} * Pk(idx_base + i);
                    end
                end
            end
            
            switch k 
                case 0
                    mpc.beq(1:mpc.nx) = -Ak*x_prev-Bdk*d(1:mpc.nd);            
                otherwise        
                    mpc.beq(k*mpc.nx+1:(k+1)*mpc.nx) = ...
                                        -Bdk*d(k*mpc.nd+1:(k+1)*mpc.nd);
            end
        end
    end
end