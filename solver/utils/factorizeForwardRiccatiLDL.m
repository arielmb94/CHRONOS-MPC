function [LDL_ric,solve_rhs_ric] = factorizeForwardRiccatiLDL( ...
    LDL_ric,solve_rhs_ric,stage,initial_stage)
coder.inline('always');

nu = size(LDL_ric,1);
last_col = size(solve_rhs_ric,2);
first_col = 1;
if initial_stage
    first_col = last_col;
end

% Factor the lower triangle; store reciprocal pivots on the diagonal.
for i = 1:nu
    pivot_reciprocal = 1 / LDL_ric(i,i,stage);
    LDL_ric(i,i,stage) = pivot_reciprocal;
    for j = i+1:nu
        Lji = LDL_ric(j,i,stage) * pivot_reciprocal;
        % Consume raw Dii*Lri before replacing this column entry with Lji.
        for r = j:nu
            LDL_ric(r,j,stage) = LDL_ric(r,j,stage) - LDL_ric(r,i,stage)*Lji;
        end
        LDL_ric(j,i,stage) = Lji;
    end
end

% Solve L*T = RHS in place, writing each nontrivial row once.
for col = first_col:last_col
    for i = 2:nu
        correction = 0;
        for j = 1:i-1
            correction = correction + LDL_ric(i,j,stage)*solve_rhs_ric(j,col);
        end
        solve_rhs_ric(i,col) = solve_rhs_ric(i,col) - correction;
    end
end
end
