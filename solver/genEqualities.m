function mpc = genEqualities(mpc,A0,B0,N,N_h_ctr,nx,nu,Pk,n_rho,varargin)
for k = 0:N-1
    % Initialize Ak and Bk as the nominal matrices for the k-th step
    Ak = A0;
    Bk = B0;
    
    % If LPV parameters are provided, compute Ak and Bk for the current step
    if nargin > 7 && ~isempty(Pk)
        idx_base = k * n_rho;
        for i = 1:n_rho
            Ak = Ak + varargin{i} * Pk(idx_base + i);
            Bk = Bk + varargin{i + n_rho} * Pk(idx_base + i);
        end
    end

    % Construct the Aeq matrix using the time-varying Ak and Bk
    if k == 0
        mpc.Aeq(1:nx,1:nu+nx) = [Bk -eye(nx)];
    
    elseif k < N_h_ctr
        mpc.Aeq(k*nx+1:(k+1)*nx,(nu+nx)*k+1-nx:(nu+nx)*(k+1)) = [Ak Bk -eye(nx)];
    else
        mpc.Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+1:(nx+nu)*(N_h_ctr-1)+nu) = ...
            Bk;
        mpc.Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr))+1:...
            (nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr-2))) = ...
            [Ak -eye(nx)];
    end
end
end