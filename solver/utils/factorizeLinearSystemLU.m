function [L, U] = factorizeLinearSystemLU(A)
% FACTORIZELINEARSYSTEMLU Compute separate Doolittle factors A = L*U.
%   A must be square and admit LU without zero pivots; no pivoting or checks.

    n = size(A, 1);
    L = eye(n);
    U = zeros(n);

    for i = 1:n
        for j = i:n
            U(i, j) = A(i, j) - L(i, 1:i-1) * U(1:i-1, j);
        end

        for j = i+1:n
            L(j, i) = (A(j, i) - L(j, 1:i-1) * U(1:i-1, i)) / U(i, i);
        end
    end
end
