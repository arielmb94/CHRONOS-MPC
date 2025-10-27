function feas = fi_box_is_feasible(cnstr)

    if any(cnstr.fi_min_x0>0) || any(cnstr.fi_max_x0>0)
        feas = 0;
    else
        feas = 1;
    end
end