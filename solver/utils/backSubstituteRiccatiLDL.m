function delta_u = backSubstituteRiccatiLDL(delta_u,LDL_ric,stage)
coder.inline('always');

nu = size(delta_u,1);
% Solve L'*delta_u = signed RHS in the destination control column.
for i = nu-1:-1:1
    correction = 0;
    for j = i+1:nu
        correction = correction + LDL_ric(j,i,stage)*delta_u(j,stage);
    end
    delta_u(i,stage) = delta_u(i,stage) - correction;
end
end
