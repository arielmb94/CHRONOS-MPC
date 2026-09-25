function Ai = appendGeneralizedConstraint(Ai, row, s_col, su_col, u_col, C_mat, Dsu_mat, D_mat)
% APPENDGENERALIZEDCONSTRAINT Inserts a unified linear constraint block.
% Inserts the Jacobian blocks for C*s + Dsu*su + D*u.

    % 1. State term (C * s_k)
    if ~isempty(s_col) && ~isempty(C_mat)
        Ai(row, s_col) = C_mat;
    end

    % 2. Rate/Past Input term (Dsu * su_k)
    if ~isempty(su_col) && ~isempty(Dsu_mat)
        Ai(row, su_col) = Dsu_mat;
    end
    
    % 3. Control Input term (D * u_k)
    if ~isempty(u_col) && ~isempty(D_mat)
        Ai(row, u_col) = D_mat;
    end
    
end
