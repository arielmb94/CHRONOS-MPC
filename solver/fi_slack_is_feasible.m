% Check feasibility of slack positivity constraints
function feas = fi_slack_is_feasible(cnstr)

    feas = 1;
    if cnstr.min_limit
        if any(cnstr.min_slack_positivity_fi_x0>0)
            feas = 0;
        end
    end
    if feas && cnstr.max_limit
        if any(cnstr.max_slack_positivity_fi_x0>0)
            feas = 0;
        end
    end

end