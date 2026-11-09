function x = solveFactoredLinearSystemLU(L, U, b)
% SOLVEFACTOREDLINEARSYSTEMLU Solve L*U*x = b for a matrix right-hand side.
%   L and U are compatible factors from FACTORIZELINEARSYSTEMLU.

    [n, m] = size(b);
    x = zeros(n, m);
    for col = 1:m
        y = forwardSubstitution(L, b(:, col));
        x(:, col) = backSubstitution(U, y);
    end
end
function y = forwardSubstitution(L, b)
    % Solves Ly = b for y (L is lower triangular).
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
    % Solves Ux = y for x (U is upper triangular).
    %   U: Upper triangular matrix
    %   y: Right-hand side vector
    %   x: Solution vector

    n = length(y);
    x = zeros(n, 1);

    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end
