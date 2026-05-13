% check feasibility of box constraints
function feas = fi_box_is_feasible(cnstr)

    feas = 1;
    if cnstr.min_limit
        if any(cnstr.fi_min_x0>0)
            feas = 0;
        end
    end
    if feas && cnstr.max_limit
        if any(cnstr.fi_max_x0>0)
            feas = 0;
        end
    end

end