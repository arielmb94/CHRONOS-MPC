function [delta_x, delta_y] = kkt_solve_combined_fast(phi, b, n, mpc_maps)
    %#codegen
    
    m = size(phi, 1) - n;
    H  = phi(1:n, 1:n);
    AT = phi(1:n, n+1:end);
    rd = b(1:n);
    rp = b(n+1:n+m);
    
    opts_LT = struct('LT', true, 'UT', false);
    opts_UT = struct('LT', false, 'UT', true);
    
    % perform a custom cholesky (using map)
    L_H = choleskyFactorization(H, mpc_maps.L_H_idx, mpc_maps.L_H_len);
    
    % replace linsolve by a custom LDL
    U = zeros(n, m);
    
    for k = 1:m
        U(:, k) = forwardSubstitutionMap(L_H, AT(:, k), mpc_maps.fwd_H_idx, mpc_maps.fwd_H_len);
    end
    
    % Schur complement
    S = U' * U;
    
    % custom cholesky
    L_S = choleskyFactorization(S, mpc_maps.L_S_idx, mpc_maps.L_S_len);
    
    % Online solve steps
    % replace linsolve by custom LDL
    u1      = forwardSubstitutionMap(L_H, rd, mpc_maps.fwd_H_idx, mpc_maps.fwd_H_len);
    by      = rp - U' * u1;
    
    % Remains the same
    u2      = linsolve(L_S,  by,      opts_LT);
    delta_y = linsolve(L_S', u2,      opts_UT);
    
    % replace by custom LDL
    u3      = -u1 - U * delta_y;
    delta_x = backSubstitutionMap(L_H', u3, mpc_maps.bwd_H_idx, mpc_maps.bwd_H_len);
end

% Custom cholesky
function L = choleskyFactorization(A, map_idx, map_len)
    n = size(A, 1);
    L = zeros(n, n);
    
    for j = 1:n
        % Diagonal
        if j == 1
            L(j,j) = sqrt(A(j,j));
        else
            L(j,j) = sqrt(A(j,j) - sum(L(j, 1:j-1).^2));
        end
        
        % number of non null elements in this column
        num_nonzeros = map_len(j);
        
        for k = 1:num_nonzeros
            i = map_idx(k, j);
            L(i, j) = (A(i, j) - L(i, 1:j-1) * L(j, 1:j-1)') / L(j, j);
        end
    end
end

function y = forwardSubstitutionMap(L, b, map_idx, map_len)
    n = length(b);
    y = zeros(n, 1);
    for i = 1:n
        num_cols = map_len(i);
        
        if num_cols == 0
            y(i) = b(i) / L(i,i);
        else
            cols = map_idx(1:num_cols, i); 
            y(i) = (b(i) - L(i, cols) * y(cols)) / L(i,i);
        end
    end
end

function x = backSubstitutionMap(U, y, map_idx, map_len)
    n = length(y);
    x = zeros(n, 1);
    for i = n:-1:1
        num_cols = map_len(i);
        
        if num_cols == 0
            x(i) = y(i) / U(i,i);
        else
            cols = map_idx(1:num_cols, i);
            x(i) = (y(i) - U(i, cols) * x(cols)) / U(i,i);
        end
    end
end