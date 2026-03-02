function feas = fi_box_is_feasible(cnstr)

    feas = 1;
    if any(cnstr.fi_min_x0>0) || any(cnstr.fi_max_x0>0)
        feas = 0;
    end

    if feas && cnstr.min_slack_nv
        if any(cnstr.min_slack_positivity_fi_x0>0)
            feas = 0;
        end
    end
    if feas && cnstr.max_slack_nv
        if any(cnstr.max_slack_positivity_fi_x0>0)
            feas = 0;
        end
    end

end