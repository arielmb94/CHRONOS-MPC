function fi_ter_x0 = get_terConst_val(x_ref,s_ter,P,v_feas)

% compute error at terminal state
x_err = x_ref-s_ter;

% Terminal Set Constraint Evaluation
fi_ter_x0 = x_err'*P*x_err-1;

% for feasibility solver: fi - v_feas <= 0
if v_feas
    fi_ter_x0 = fi_ter_x0 - v_feas;
end

end