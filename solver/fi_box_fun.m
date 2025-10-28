function cnstr = fi_box_fun(cnstr,x,N,n,v_feas)

if ~isempty(cnstr.min)

    % If active-set optimization enabled reset indicator vector
    if cnstr.min_activ_set
        cnstr.min_activ_indicator(:) = 0;
    end

    for k = 1:N/n

        cnstr.fi_min_x0((k-1)*n+1:k*n) = cnstr.min-x((k-1)*n+1:k*n);

        % check wich constraints are active
        if cnstr.min_activ_set

            % check wich conditions enabled for the k horizon step
            step_indicator = x((k-1)*n+1:k*n) < cnstr.min_activ_lim;
            % update the active condition indicator vector
            cnstr.min_activ_indicator = cnstr.min_activ_indicator | step_indicator;

        end
    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        cnstr.fi_min_x0 = cnstr.fi_min_x0 - v_feas;
    end

end

if ~isempty(cnstr.max)

    % If active-set optimization enabled reset indicator vector
    if cnstr.max_activ_set
        cnstr.max_activ_indicator(:) = 0;
    end

    for k = 1:N/n

        cnstr.fi_max_x0((k-1)*n+1:k*n) = x((k-1)*n+1:k*n)-cnstr.max;

        % check wich constraints are active
        if cnstr.max_activ_set

            % check wich conditions enabled for the k horizon step
            step_indicator = x((k-1)*n+1:k*n) > cnstr.max_activ_lim;
            % update the active condition indicator vector
            cnstr.max_activ_indicator = cnstr.max_activ_indicator | step_indicator;

        end
    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        cnstr.fi_max_x0 = cnstr.fi_max_x0 - v_feas;
    end

end

end