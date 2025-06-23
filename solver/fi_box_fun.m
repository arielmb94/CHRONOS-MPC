function [fi_low_x0,fi_up_x0] = fi_box_fun(x,low_lim,up_lim,N,n,v_feas)

if ~isempty(low_lim)

    fi_low_x0 = zeros(N,1);
    for k = 1:N/n

        fi_low_x0((k-1)*n+1:k*n) = low_lim-x((k-1)*n+1:k*n);

    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        fi_low_x0 = fi_low_x0 - v_feas;
    end

else

    fi_low_x0 = [];

end

if ~isempty(up_lim)

    fi_up_x0 = zeros(N,1);
    for k = 1:N/n

        fi_up_x0((k-1)*n+1:k*n) = x((k-1)*n+1:k*n)-up_lim;

    end

    % for feasibility solver: fi - v_feas <= 0
    if v_feas
        fi_up_x0 = fi_up_x0 - v_feas;
    end

else

    fi_up_x0 = [];

end

end