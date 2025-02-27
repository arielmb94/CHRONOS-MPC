function x = solveLinearSystemLU(A, b)
    % "solveLinearSystemLU" Solves Ax = b using LU decomposition.
    %   A: Coefficient matrix
    %   b: Right-hand side vector
    %   x: Solution vector
    
    [L, U] = luDecomposition(A);
    y = forwardSubstitution(L, b);
    x = backSubstitution(U, y);
end
function [L, U] = luDecomposition(A)
    % "luDecomposition" Performs LU decomposition without pivoting.
    %   A: Coefficient matrix
    %   L: Lower triangular matrix
    %   U: Upper triangular matrix
    
    n = size(A, 1);
    L = eye(n);
    U = zeros(n);
    
    for i = 1:n
        % Upper triangular matrix U
        for j = i:n
            U(i, j) = A(i, j) - L(i, 1:i-1) * U(1:i-1, j);
        end
        
        % Lower triangular matrix L
        for j = i+1:n
            L(j, i) = (A(j, i) - L(j, 1:i-1) * U(1:i-1, i)) / U(i, i);
        end
    end
end
function y = forwardSubstitution(L, b)
    % "forwardSubstitution" Solves Ly = b for y (L is lower triangular).
    %   L: Lower triangular matrix
    %   b: Right-hand side vector
    %   y: Solution vector
    
    n = length(b);
    y = zeros(n, 1);
    
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1));
    end
end
function x = backSubstitution(U, y)
    % "backSubstitution" Solves Ux = y for x (U is upper triangular).
    %   U: Upper triangular matrix
    %   y: Right-hand side vector
    %   x: Solution vector
    
    n = length(y);
    x = zeros(n, 1);
    
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end

