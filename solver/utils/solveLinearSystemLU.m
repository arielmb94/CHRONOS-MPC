function x = solveLinearSystemLU(A, b)
% SOLVELINEARSYSTEMLU Factor A and solve A*x = b without pivoting.

    [L, U] = factorizeLinearSystemLU(A);
    x = solveFactoredLinearSystemLU(L, U, b);
end
