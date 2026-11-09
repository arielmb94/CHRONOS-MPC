function Ai = appendGeneralizedConstraint(Ai, row, x_col, xu_col, u_col, C_mat, Dxu_mat, D_mat)
% APPENDGENERALIZEDCONSTRAINT Inserts a unified linear constraint block.
% Inserts the Jacobian blocks for C*x + Dxu*xu + D*u.

    % 1. State term (C * x_k)
    if ~isempty(x_col) && ~isempty(C_mat)
        Ai(row, x_col) = C_mat;
    end

    % 2. Rate/Past Input term (Dxu * xu_k)
    if ~isempty(xu_col) && ~isempty(Dxu_mat)
        Ai(row, xu_col) = Dxu_mat;
    end
    
    % 3. Control Input term (D * u_k)
    if ~isempty(u_col) && ~isempty(D_mat)
        Ai(row, u_col) = D_mat;
    end
    
end
