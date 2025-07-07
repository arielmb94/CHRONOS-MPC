function cnstr = fi_box_fun(cnstr,x,N,n,v_feas)

if ~isempty(cnstr.min)

    for k = 1:N/n

        cnstr.fi_min_x0((k-1)*n+1:k*n) = cnstr.min-x((k-1)*n+1:k*n);

    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        cnstr.fi_min_x0 = cnstr.fi_min_x0 - v_feas;
    end

end

if ~isempty(cnstr.max)

    for k = 1:N/n

        cnstr.fi_max_x0((k-1)*n+1:k*n) = x((k-1)*n+1:k*n)-cnstr.max;

    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        cnstr.fi_max_x0 = cnstr.fi_max_x0 - v_feas;
    end

end

end