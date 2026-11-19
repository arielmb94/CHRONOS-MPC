function scaled_rhs_ric = scaleRiccatiLDLRHS( ...
    scaled_rhs_ric,solve_rhs_ric,LDL_ric,k)
coder.inline('always');

% Retain Z = D^(-1)*T for the Schur correction and forward recovery.
for col = 1:size(solve_rhs_ric,2)
    for i = 1:size(solve_rhs_ric,1)
        scaled_rhs_ric(i,col,k) = LDL_ric(i,i,k+1)*solve_rhs_ric(i,col);
    end
end
end
