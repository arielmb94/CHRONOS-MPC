% check feasibility of box constraints
function feas = fi_box_is_feasible(cnstr)

    feas = true;
    if cnstr.min_limit
        if any(cnstr.fi_min_x0>0)
            feas = false;
        end
    end
    if feas && cnstr.max_limit
        if any(cnstr.fi_max_x0>0)
            feas = false;
        end
    end

end