function X = solveLinearSystemLDL(A, B)
    % Solves AX = B using LDL' factorization for symmetric matrices.
    %   A: Symmetric matrix (n x n)
    %   B: Right-hand side vector or matrix (n x m)
    %   X: Solution vector or matrix (n x m)

    [L, d] = ldlFactorization(A);
    Y = forwardSubstitution(L, B);
    Z = diagonalSolve(d, Y);
    X = backSubstitutionTranspose(L, Z);
end

function [L, d] = ldlFactorization(A)
    n = size(A, 1);
    L = eye(n);
    d = zeros(n, 1);
    
    for i = 1:n
        % Compute D(i) = A(i,i) - sum_{k=1}^{i-1} L(i,k)^2 * d(k)
        sum_d = 0;
        for k = 1:i-1
            sum_d = sum_d + L(i, k) * L(i, k) * d(k);
        end
        d(i) = A(i, i) - sum_d;
        
        for j = i+1:n
            sum_l = 0;
            for k = 1:i-1
                sum_l = sum_l + L(j, k) * L(i, k) * d(k);
            end
            L(j, i) = (A(j, i) - sum_l) / d(i);
        end
    end
end

function Y = forwardSubstitution(L, B)
    % Solves L * Y = B for Y (L is unit lower triangular)
    [n, m] = size(B);
    Y = zeros(n, m);
    
    for col = 1:m
        for i = 1:n
            sum_y = 0;
            for k = 1:i-1
                sum_y = sum_y + L(i, k) * Y(k, col);
            end
            Y(i, col) = B(i, col) - sum_y;
        end
    end
end

function Z = diagonalSolve(d, Y)
    % Solves diag(d) * Z = Y
    [n, m] = size(Y);
    Z = zeros(n, m);
    
    for col = 1:m
        for i = 1:n
            Z(i, col) = Y(i, col) / d(i);
        end
    end
end

function X = backSubstitutionTranspose(L, Z)
    % Solves L' * X = Z for X without explicitly transposing L in memory
    [n, m] = size(Z);
    X = zeros(n, m);
    
    for col = 1:m
        for i = n:-1:1
            sum_x = 0;
            for k = i+1:n
                % (L')_{i,k} corresponds to L(k, i)
                sum_x = sum_x + L(k, i) * X(k, col);
            end
            X(i, col) = Z(i, col) - sum_x;
        end
    end
end